#!/bin/bash
# GSS bitset micro-PR benchmark kit: main vs #135 (iteration) vs #136 (Hall early stop) vs
# #137 (recount only filtered domains) vs three fused-popcount variants, each built with GCC
# and clang (Release, so -O3 -march=native). See README.md and issue #138.
#
# usage: ./run-kit.bash [REPS]        (default 4)
# env:   CORE=n      core to pin to on Linux (default 2; on a hybrid Intel CPU, pick a P-core)
#        CXX_GCC / CXX_CLANG   compilers (default g++ / clang++; on macOS g++ is usually clang,
#                              so set CXX_GCC=g++-15 or similar)
#        JOBS=n      parallel build jobs (default: all cores)
#        ARMS="..."  subset of arms to build/run (default: all fourteen)
#        SMOKE=1     quick check: three short instances, one rep, separate results file
#
# Runs one solver at a time, single-threaded, pinned where possible. Writes
# results-<host>.txt: arm cfg instance rep padlen instructions cycles runtime_ms nodes status rc
# and machine-<host>.txt describing the CPU and compilers. It is resumable: rerun the same
# command after an interruption and it skips runs already in the results file.
set -u
here=$(cd "$(dirname "$0")" && pwd)
reps=${1:-4}
host=$(hostname -s)
work=$here/work; mkdir -p "$work"
jobs=${JOBS:-$( (nproc || sysctl -n hw.ncpu) 2>/dev/null)}
repo=https://github.com/ciaranm/glasgow-subgraph-solver.git
os=$(uname -s)

# The exact commits measured. Pinned, so that later pushes to the PR branches or to main do
# not change what a run measures; every machine builds the same code.
sha_main=01ef6caec4dbf0ef0206689d769f0b67aacd38b5   # main
sha_pr1=937a632c48b8f1b6ea15c7b16504ac12100b1ab7    # #135 svo-bitset-iteration
sha_pr2=bddb4e581d8127769038153daae03ad44a7a0e41    # #136 cheap-all-different-hall-early-stop
sha_pr3=71108cc9d7ed28e76508f84c7d01015d3310e8f8    # #137 recount-only-filtered-domains

suffix=""; instances_file=$here/instances.txt
if [ "${SMOKE:-0}" = 1 ]; then
    suffix=-smoke; reps=1
    instances_file=$work/smoke-instances.txt
    grep -E ' (si6_b06_m800\.02|all-images-CVIU11-pattern12-target88|all-images-CVIU11-pattern24-target22) ' \
        "$here/instances.txt" > "$instances_file"
fi
out=$here/results-$host$suffix.txt

have() { command -v "$1" > /dev/null 2>&1; }
prefix=(); perf_ok=0
if [ "$os" = Linux ]; then
    have taskset && prefix+=(taskset -c "${CORE:-2}")
    have setarch && prefix+=(setarch "$(uname -m)" -R)
    have perf && perf stat -x, -e instructions:u,cycles:u -o /dev/null true 2> /dev/null && perf_ok=1
fi

# --- machine description (no hostname: this gets posted publicly) -----------------------
{
    uname -srm
    if [ "$os" = Darwin ]; then
        sysctl -n machdep.cpu.brand_string
        sysctl hw.perflevel0.physicalcpu hw.perflevel1.physicalcpu hw.memsize 2>/dev/null
        sw_vers 2>/dev/null
        pmset -g batt 2>/dev/null | head -1
    else
        lscpu
        cat /sys/devices/system/cpu/cpufreq/boost 2>/dev/null | sed 's/^/boost /'
        cat /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null | sed 's/^/intel_pstate no_turbo /'
        cat /sys/devices/system/cpu/cpu${CORE:-2}/cpufreq/scaling_governor 2>/dev/null | sed 's/^/governor /'
        cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null | sed 's/^/perf_event_paranoid /'
        grep -m1 '^microcode' /proc/cpuinfo 2>/dev/null
        echo "pinned core ${CORE:-2}: $(have taskset && echo taskset || echo 'not pinned (no taskset)')"
    fi
    echo "perf counters: $( ((perf_ok)) && echo yes || echo 'no (runtime only)')"
    echo "gcc arm compiler ${CXX_GCC:-g++}: $(${CXX_GCC:-g++} --version 2>&1 | head -1)"
    echo "clang arm compiler ${CXX_CLANG:-clang++}: $(${CXX_CLANG:-clang++} --version 2>&1 | head -1)"
    echo "cmake: $(cmake --version 2>&1 | head -1)"
    echo "commits: main $sha_main pr1 $sha_pr1 pr2 $sha_pr2 pr3 $sha_pr3"
    echo "kit: $(git -C "$here" rev-parse HEAD 2>/dev/null || echo unknown)"
} > "$here/machine-$host$suffix.txt" 2>&1

# --- sources ---------------------------------------------------------------------------
if [ ! -d "$work/src" ]; then
    git clone -q "$repo" "$work/src" || exit 1
fi
srcdir() { echo "$work/src-$1"; }
for a in main pr1 pr2 pr3 v; do
    case $a in
        pr1) r=$sha_pr1 ;;
        pr2) r=$sha_pr2 ;;
        pr3) r=$sha_pr3 ;;
        *) r=$sha_main ;;
    esac
    d=$(srcdir $a)
    [ -d "$d" ] && continue
    git -C "$work/src" cat-file -e "$r^{commit}" 2> /dev/null || git -C "$work/src" fetch -q origin "$r" || exit 1
    git -C "$work/src" worktree add -q --detach "$d" "$r" || exit 1
    if [ $a = v ]; then git -C "$d" apply "$here/fused-variants.patch" || exit 1; fi
done

# --- builds ----------------------------------------------------------------------------
all_arms="main-g main-c pr1-g pr1-c pr2-g pr2-c pr3-g pr3-c v1-g v1-c v2-g v2-c v3-g v3-c"
arms=${ARMS:-$all_arms}
mkdir -p "$work/bin"
for arm in $arms; do
    [ -x "$work/bin/$arm" ] && continue
    base=${arm%-*}; comp=${arm##*-}
    cxx=${CXX_GCC:-g++}; [ $comp = c ] && cxx=${CXX_CLANG:-clang++}
    flags=""; s=$(srcdir $base)
    case $base in v1|v2|v3) flags="-DGSS_COUNT_VARIANT=${base#v}"; s=$(srcdir v) ;; esac
    b=$work/build-$arm
    echo "building $arm with $cxx" >&2
    cmake -S "$s" -B "$b" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$cxx" -DCMAKE_CXX_FLAGS="$flags" > "$b.log" 2>&1 &&
        cmake --build "$b" -j "$jobs" --target glasgow_subgraph_solver >> "$b.log" 2>&1 ||
        { echo "build of $arm failed, see $b.log"; exit 1; }
    cp "$b/glasgow_subgraph_solver" "$work/bin/$arm"
done

# --- runner ----------------------------------------------------------------------------
touch "$out"
inst=$here/instances
run_one() { # arm cfg name pat tgt rep
    local arm=$1 cfg=$2 name=$3 pat=$4 tgt=$5 rep=$6 tmp args pad ins=NA cyc=NA
    grep -q "^$arm $cfg $name $rep " "$out" && return   # resumable
    args=(--format lad); [ $cfg = nosupp ] && args+=(--no-supplementals)
    # The environment's size shifts the initial stack; vary it by rep so that one unlucky
    # stack alignment does not bias every rep of an arm the same way.
    pad=$(printf '%*s' $(( (rep % 4) * 97 + 1 )) '' | tr ' ' x)
    tmp=$(mktemp -d)
    if ((perf_ok)); then
        env -i PATH=/usr/bin:/bin PAD=$pad GLIBC_TUNABLES=glibc.malloc.mmap_threshold=33554432:glibc.malloc.trim_threshold=4294967295 \
            "${prefix[@]}" perf stat -x, -e instructions:u,cycles:u -o $tmp/perf \
            "$work/bin/$arm" "${args[@]}" --timeout 300 "$inst/$pat" "$inst/$tgt" > $tmp/out 2>&1
    else
        env -i PATH=/usr/bin:/bin PAD=$pad GLIBC_TUNABLES=glibc.malloc.mmap_threshold=33554432:glibc.malloc.trim_threshold=4294967295 \
            "${prefix[@]}" "$work/bin/$arm" "${args[@]}" --timeout 300 "$inst/$pat" "$inst/$tgt" > $tmp/out 2>&1
    fi
    local rc=$?
    if ((perf_ok)); then ins=$(awk -F, '/instructions/ {print $1}' $tmp/perf); cyc=$(awk -F, '/cycles/ {print $1}' $tmp/perf); fi
    local rt=$(awk '/^runtime = / {print $3}' $tmp/out) nodes=$(awk '/^nodes = / {print $3}' $tmp/out) st=$(awk '/^status = / {print $3}' $tmp/out)
    echo "$arm $cfg $name $rep ${#pad} ${ins:-NA} ${cyc:-NA} ${rt:-NA} ${nodes:-NA} ${st:-NA} $rc" >> "$out"
    rm -rf $tmp
}
# Interleave arms within each (instance, rep), rotating their order, so slow drift
# (thermals, background load) hits every arm alike.
total=$(( $(wc -l < "$instances_file") * reps ))
k=0
for rep in $(seq 0 $((reps - 1))); do
    while read -r cfg name pat tgt; do
        k=$((k + 1)); printf '\r%d/%d' $k $total >&2
        for arm in $( printf '%s\n' $arms | awk -v s=$((k + rep)) '{a[NR]=$0} END {for (i = 0; i < NR; ++i) print a[(i + s) % NR + 1]}' ); do
            run_one $arm $cfg $name $pat $tgt $rep
        done
    done < "$instances_file"
done
echo >&2; echo "done: $out"
