# this script is here to make my life easier when testing multiple instances
struct Options
    ins::String     # one instance name (leave blank if you want all instances)
    insid::String   # instance id
    bench::String   # directory containing instances.
    pbopath::String   # directory containing veripb.
    solveurpath::String   # directory containing solver.
    proofs::String   # directory containing solver.
    veripb::Bool    # veripb comparison (need veripb3.0 installed)
    trace::Bool     # veripb trace
    profiling::Bool # profiling
    timelimit::Int # time limit
    rand::Bool # randomize the instances
    pattern::String # pattern name
    target::String  # target name
    format::String  # format of the instances ex lad or directedlad
end
function parseargs(args)
    ins = ""
    # bench = "/home/arthur_gla/veriPB/newSIPbenchmarks"
    bench = "/users/grad/arthur/newSIPbenchmarks"
    # pbopath = "/home/arthur_gla/veriPB/subgraphsolver/veripb-dev"
    pbopath = "/users/grad/arthur/pboxide-dev"
    # solveurpath = "/home/arthur_gla/veriPB/subgraphsolver/glasgow-subgraph-solver/build"    
    solveurpath = "/users/grad/arthur/glasgow-subgraph-solver/build"
    # proofs = "/home/arthur_gla/veriPB/subgraphsolver/proofs"
    # proofs = "/home/arthur_gla/veriPB/subgraphsolver/nolabelsproofs3"
    # proofs = "/scratch/arthur/proofs_test_veriPBtrim_mem_out"
    # proofs = "users/grad/arthur/proofs"
    proofs = "users/grad/arthur/proofs2"
    veripb = false
    trace = false
    prof = false
    rand = false
    format = "lad"
    insid = pattern = target = ""
    tl = 600
    for (i, arg) in enumerate(args)
        if arg == "cd" cd() end # hack to add cd in paths
        if arg in ["--trace","-trace","trace"] trace = true end
        if arg in ["--profiling","-profiling","profiling","prof","-prof","--prof"] prof = true end
        if arg in ["--trace","-trace","trace","-tr","tr"] trace = true end
        if arg in ["rand","random"] rand = true end
        if arg in ["veripb","verif"] veripb = true end
        if arg in ["timelimit","tl"] tl = parse(Int, args[i+1]) end
        if arg in ["insid","ins"] insid = args[i+1] end
        if arg in ["pattern","p"] pattern = args[i+1] end
        if arg in ["target","t"] target = args[i+1] end
        if arg in ["format"] format = args[i+1] end
        if arg in ["directed"] format = "directedlad" end
    end
    return Options(ins,insid,bench,pbopath,solveurpath,proofs,veripb,trace,prof,tl,rand,pattern,target,format)
end
const CONFIG = parseargs(ARGS)
const benchs = CONFIG.bench 
const pbopath = CONFIG.pbopath 
const solver = CONFIG.solveurpath 
const proofs = CONFIG.proofs 
const tl = CONFIG.timelimit
const format = CONFIG.format
const extention = ".pbp"
function main()
    cd(solver)
    try run(`make`)
    catch e println(e) end
    if "bio" in ARGS
    run_bio_solver()
    else 
    run_LV_solver()
    end
end
function run_bio_solver()
    path = string(benchs,"/biochemicalReactions")
    cd()
    graphs = cd(readdir, path)
    n = length(graphs)
    if CONFIG.insid != "" || CONFIG.pattern != "" && CONFIG.target != ""
        pattern = length(CONFIG.pattern)>0 && CONFIG.pattern[end-2:end]=="txt" ? CONFIG.pattern[1:end-4] : CONFIG.pattern
        target = length(CONFIG.target)>0 && CONFIG.target[end-2:end]=="txt" ? CONFIG.target[1:end-4] : CONFIG.target
        ins = CONFIG.insid
        if ins != ""
            ins = ins[1:3]=="bio" ? ins : string("bio",ins)
            pattern = CONFIG.insid[4:6]
            target = CONFIG.insid[7:9]
        else
            ins = string("bio",pattern,target)
        end
        solve(ins,path,pattern*".txt",path,target*".txt",0,200_000_000_000)
        if CONFIG.veripb
            cd()
            cd(CONFIG.pbopath)
            runpboxide(ins)
        end
        println()
    else
    t=0
    stats = [stat(path*'/'*file).size for file in graphs]
    println("stats: ", stats)
    p = sortperm(stats)
    for i in eachindex(graphs), j in eachindex(graphs)
        if CONFIG.rand
            i = rand(1:length(graphs))
            j = rand(1:length(graphs))
        end
        pattern = graphs[p[i]]
        target = graphs[p[j]]
        # t+=1
        # if t>30 break end
        theokbigs = ["bio055013","bio064013","bio066013","bio122013","bio123013","bio002014"]

        if pattern != target

            ins = string("bio",pattern[1:end-4],target[1:end-4])
            if ins in ["bio007013"] continue end # skip this instance, it is too big
            solve(ins,path,pattern,path,target)
            if isfile(String("$proofs/$ins$extention"))
                res = read(`tail -n 2 $proofs/$ins$extention`,String)
                if length(res)<16 || res[1:16] != "conclusion UNSAT"
                    printstyled("sat or unfinished proof\n", color=:green)
                else
                    if CONFIG.veripb
                        cd()
                        cd(CONFIG.pbopath)
                        runpboxide(ins)
                    end
                end
            end
            println()
        end
    end end
end
function run_LV_solver()
    path = string(benchs,"/LV")
    cd()
    graphs = cd(readdir, path)
    n = length(graphs)
    if CONFIG.pattern != "" && CONFIG.target != "" || CONFIG.insid != ""
        ins = string("LV",CONFIG.insid)
        i = findlast(c -> c == 'g', CONFIG.insid)
        pattern = CONFIG.insid[1:i-1]
        target = CONFIG.insid[i:end]
        solve(ins,path,pattern,path,target,0,200_000_000_000)
        if CONFIG.veripb
            cd()
            cd(CONFIG.pbopath)
            runpboxide(ins)
        end
        println()
    else
    t=0
    stats = [stat(path*'/'*file).size for file in graphs]
    println("stats: ", stats)
    p = sortperm(stats)
    for i in eachindex(graphs), j in eachindex(graphs)
        if CONFIG.rand
            i = rand(1:length(graphs))
            j = rand(1:length(graphs))
        end
        pattern = graphs[p[i]]
        target = graphs[p[j]]
    # for target in graphs[1:end], pattern in graphs[1:end]
        if pattern != target && !(pattern in ["g2","g3"])
            ins = string("LV",pattern,target)
            solve(ins,path,pattern,path,target)
            if isfile(String("$proofs/$ins$extention"))
                res = read(`tail -n 2 $proofs/$ins$extention`,String)
                if length(res)<16 || res[1:16] != "conclusion UNSAT"
                    printstyled("sat or unfinished proof\n", color=:green)
                else
                    if CONFIG.veripb
                        cd()
                        cd(CONFIG.pbopath)
                        runpboxide(ins)
                    end
                end
            end
            println()
        end
    end end
end
function solve(ins,pathpat,pattern,pathtar,target,minsize=2_000_000_000,maxsize=20_000_000_000,remake=true,verbose=false)
    if remake || !isfile(string(proofs,"/",ins,".opb")) || !isfile(string(proofs,"/",ins,extention)) || 
            length(read(`tail -n 1 $proofs/$ins$extention`,String)) < 24 ||
            read(`tail -n 1 $proofs/$ins$extention`,String)[1:24] != "end pseudo-Boolean proof"
        print(ins,' ')
        cd(CONFIG.solveurpath)
        # println("timeout $tl ./glasgow_subgraph_solver --prove $proofs/$ins --no-clique-detection --format $format $pathpat/$pattern $pathtar/$target")
        t = @elapsed begin
            # p = run(pipeline(`timeout $timeout ./$solver --prove $proofs/$ins --no-supplementals --no-clique-detection --format $format $pathpat/$pattern $pathtar/$target`, devnull),wait=false); wait(p)
            try
                if verbose
                    p = run(`timeout $tl ./glasgow_subgraph_solver --prove $proofs/$ins --no-clique-detection --format $format $pathpat/$pattern $pathtar/$target`)
                else
                    # println(`timeout $tl ./glasgow_subgraph_solver --prove $proofs/$ins --no-clique-detection --format $format $pathpat/$pattern $pathtar/$target`)
                    redirect_stdio(stdout = devnull) do
                    p = run(`timeout $tl ./glasgow_subgraph_solver --prove $proofs/$ins --no-clique-detection --format $format $pathpat/$pattern $pathtar/$target`)
                    # p = run(`timeout $tl ./glasgow_subgraph_solver --no-clique-detection $pathpat/$pattern $pathtar/$target`) # no proof version for comparison
                    end
                end
            catch e
            end
        end
        t+=0.01
        ok = false
        print(prettytime(t))
        size = stat(string(proofs,"/",ins,extention)).size
        if t>tl
            printstyled(" timeout "; color = :red)
        elseif read(`tail -n 2 $proofs/$ins$extention`,String)[1:14] == "conclusion SAT"
            printstyled(" sat     "; color = 166)
        elseif minsize > size            
            printstyled(" toosmal ",prettybytes(size); color = :yellow)            
        elseif maxsize < size            
            printstyled(" toobig  ",prettybytes(size); color = :red)
        else printstyled(" OK      ",prettybytes(size); color = :green)
            ok = true
            # g = ladtograph(pathpat,pattern)
            # draw(PNG(string(proofs,"/aimg/graphs/",ins,pattern[1:3],".png"), 16cm, 16cm), gplot(g))
            # g = ladtograph(pathtar,target)
            # draw(PNG(string(proofs,"/aimg/graphs/",ins,target[1:3],".png"), 16cm, 16cm), gplot(g))
        end
        if !ok
            run(`rm -f $proofs/$ins$extention`)
            run(`rm -f $proofs/$ins.opb`)
        end
    end
end
function runpboxide(file)
    t1 = t2 = 0
    prof = CONFIG.profiling ? "flamegraph" : "r"
    tll = 10*tl
    t1 = @elapsed begin
        try
            printstyled("    veriPB check ", color=:blue)
                if CONFIG.trace
                    v1 = run(`timeout $tll cargo $prof -- $proofs/$file.opb $proofs/$file$extention `)
                else
                    redirect_stdio(stdout = devnull,stderr = devnull) do
                        v1 = run(`timeout $tll cargo $prof -- $proofs/$file.opb $proofs/$file$extention `)
                    end
                end
        catch e
            printstyled("fail ", color=:red)
            # print(e)
        end
    end
    if t1>tll
        printstyled("timeout "; color = :red)
    end
    printstyled(prettytime(t1); color = :cyan)
    # printstyled(prettytime(t1),"  ",prettytime(t2),"  "; color = :cyan)
end
function prettytime(b)
    if b<0.01
        return  string(0)
    elseif b<0.1
        return  string(round(b; sigdigits=1))
    elseif b<1
        return  string(round(b; sigdigits=2))
    else
        return  string(round(b; sigdigits=3))
    end
end
function prettybytes(b)
    if b>=10^9
        return string(round(b/(10^9); sigdigits=4)," GB")
    elseif b>=10^6
        return string(round(b/(10^6); sigdigits=4)," MB")
    elseif b>=10^3
        return string(round(b/(10^3); sigdigits=4)," KB")
    else
        return  string(round(b; sigdigits=4)," B")
    end
end 
main()



#=
julia GlasgowPB3trimnalyser.jl LVg6g12 cshow
julia GlasgowPB3trimnalyser.jl LVg7g71 cshow
julia GlasgowPB3trimnalyser.jl LVg11g72 cshow
julia GlasgowPB3trimnalyser.jl LVg17g19 cshow ?
julia GlasgowPB3trimnalyser.jl LVg18g59 cshow
julia GlasgowPB3trimnalyser.jl LVg16g58 cshow
julia GlasgowPB3trimnalyser.jl LVg12g62 cshow

 julia script.jl ins g6g12 nv

julia GlasgowPB3trimnalyser.jl LVg26g100 cshow
solve 14s
LVg26g100 & 63.85 MB &          & 4.682 MB & 0    & 0    & 40.6 ( 4.97 34.4 1.22 0 ) \\\hline

LVg26g100 & 62.83 MB &          & 5.073 MB & 0    & 0    & 46.5 (5.59 39.9 0.93 0   ) \\\hline

conflict analysis
                                  4.085 MB                 48.6  5.64 42.1 0.84 0
ca + wtrim
                                  4.085 MB                 46.6  4.59 41.1 0.91 0
ca + wtrim + @labels
                                  4.982 MB                 44.5  4.37 39.2 0.93 0
remarques, inj1 n'est pas utilise mais deg36=1 l'est


berhan watch list a pleins d'elements ce qui est bizare.
mais ce qui tue le truc c'est qu'il faut parcourir toute la liste pour trouver le dernier element et le supprimer.
ca peut venir du fait que on regarde le defered set et qu'on ne compute pas le slack de ces contraintes
quand on regarde le code de veripb meme quand on suprime la database complete, on a pas ca parce que les watched literaux on ete update par des rup call.

dans la propagation tu ajoute dans les watched literal tous les lits et du coup tu galere a  les suprimer.


nice instances 
LVg3g52 & 10.92 MB & 582.9 KB & 194.6 KB & 3.25 & 4.82 & 5.67 (0.8  4.61 0.01 0.25) \\\hline



weird instances.
LVg3g17 & 1.665 MB & 499.4 KB & 643.8 KB & 0.46 & 0.79 & 1.19 (0.21 0.66 0.02 0.3 ) \\\hline
LVg2g80 & 5.219 MB & 6.455 MB & 414.1 KB & 44.2 & 12.5 & 100.0(2.7  96.5 0.16 0.86) \\\hline
LVg3g41 & 8.664 MB & 6.284 MB & 421.1 KB & 45.3 & 46.4 & 261.0(4.72 245.01.02 9.94) \\\hline

LVg3g39 & 6.552 MB & 4.058 MB & 328.6 KB & 30.0 & 28.4 & 168.0(3.27 157.00.88 6.79) \\\hline

bio027085 & 2.922 MB & 333.2 KB & 75.15 KB & 4.28 & 32.4 & 3.7  (2.35 0.98 0.02 0.35) \\\hline

bio068087 & 6.905 MB & 6.235 MB & 292.0 KB & 9.64 & 125.0& 18.3 (8.09 9.37 0.05 0.81) \\\hline
bio084056 & 12.74 MB &          &          & 21.9 &      &      (                   ) \\\hline


bio169086 & 16.12 MB & 7.491 MB & 106.5 KB & 16.4 & 126.0& 23.6 (11.2 10.0 1.79 0.57) \\\hline

bio092151 & 14.28 MB & 22.42 MB & 1.194 MB & 18.8 & 871.0& 73.2 (23.2 47.9 0.11 1.99) \\\hline

bio171015 & 14.93 MB & 27.89 MB & 1.092 MB & 146.0& 268.0& 1500 (9.35 1370.8.95 114) \\\hline
bio170075 & 18.54 MB & 14.15 MB & 585.5 KB & 192.0& 173.0& 571.0(15.3 523.02.26 30.8) \\\hline

bio068151 & 19.0 MB  &          &          & 44.1^C

a surveiller bio037111 elle parais hyper longue.



bio029051 & 18.88 MB &          & 940.9 KB & 36.5 & 0    & 82.4 (33.3 46.9 0.2  2.04) \\\hline
bio170075 & 19.07 MB &          & 609.0 KB & 158.0& 0    & 564.0(14.6 510.02.2  37.5) \\\hline

interessante ?
bio116023


LVg46g74 10.1 OK      5.16 MB
LVg22g79 54.5 toobig  2.18 GB
LVg22g30 24.9 toobig  1.143 GB

exemple de degre non trime jusqau bout 
bio105014

./glasgow_subgraph_solver --prove /home/arthur_gla/veriPB/subgraphsolver/proofs/bio170075 --no-clique-detection --format directedlad /home/arthur_gla/veriPB/newSIPbenchmarks/biochemicalReactions/170.txt /home/arthur_gla/veriPB/newSIPbenchmarks/biochemicalReactions/075.txt
./glasgow_subgraph_solver --prove /home/arthur_gla/veriPB/subgraphsolver/proofs/LVg6g12 --no-clique-detection --format directedlad /home/arthur_gla/veriPB/newSIPbenchmarks/LV/g6 /home/arthur_gla/veriPB/newSIPbenchmarks/LV/g12
=#