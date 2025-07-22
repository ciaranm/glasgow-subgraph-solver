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
end
function parseargs(args)
    ins = ""
    bench = "/home/arthur_gla/veriPB/newSIPbenchmarks"
    pbopath = "/home/arthur_gla/veriPB/subgraphsolver/pboxide-dev"
    solveurpath = "/home/arthur_gla/veriPB/subgraphsolver/glasgow-subgraph-solver/build"
    proofs = "/home/arthur_gla/veriPB/subgraphsolver/proofs"
    veripb = true
    trace = false
    prof = false
    insid = ""
    tl = 600
    for (i, arg) in enumerate(args)
        if arg == "cd" cd() end # hack to add cd in paths
        if arg in ["--trace","-trace","trace"] trace = true end
        if arg in ["--profiling","-profiling","profiling","prof","-prof","--prof"] prof = true end
        if arg in ["--trace","-trace","trace","-tr","tr"] trace = true end
        if arg in ["noveripb","nv"] veripb = false end
        if arg in ["timelimit","tl"] tl = parse(Int, args[i+1]) end
        if arg in ["insid","ins"] insid = args[i+1] end
    end
    return Options(ins,insid,bench,pbopath,solveurpath,proofs,veripb,trace,prof,tl)
end
const CONFIG = parseargs(ARGS)
const benchs = CONFIG.bench 
const pbopath = CONFIG.pbopath 
const solver = CONFIG.solveurpath 
const proofs = CONFIG.proofs 
const tl = CONFIG.timelimit
const extention = ".pbp"
function main()
    cd(solver)
    try run(`make`)
    catch e println(e) end
    # run_bio_solver()
    run_LV_solver()
end
function run_bio_solver()
    path = string(benchs,"/biochemicalReactions")
    cd()
    graphs = cd(readdir, path)
    n = length(graphs)
    if CONFIG.insid != ""
        ins = string("bio",CONFIG.insid)
        solve(ins,path,CONFIG.insid[1:3]*".txt",path,CONFIG.insid[4:6]*".txt","lad",0,1_000_000_000)
        if CONFIG.veripb
            cd()
            cd(CONFIG.pbopath)
            runpboxide(ins)
        end
    else
    t=0
    for target in graphs[1:end], pattern in graphs[1:end]
        # target = graphs[rand(1:n)]
        # pattern = graphs[rand(1:n)]
        # t+=1
        # if t>30 break end
        theokbigs = ["bio055013","bio064013","bio066013","bio122013","bio123013","bio002014"]

        if pattern != target

            ins = string("bio",pattern[1:end-4],target[1:end-4])
            if ins in ["bio007013"] continue end # skip this instance, it is too big
            solve(ins,path,pattern,path,target,"lad")
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
        end
    end end
end
#=
julia GlasgowPB3trimnalyser.jl LVg6g12 cshow
julia GlasgowPB3trimnalyser.jl LVg7g71 cshow
julia GlasgowPB3trimnalyser.jl LVg11g72 cshow
julia GlasgowPB3trimnalyser.jl LVg17g19 cshow ?
julia GlasgowPB3trimnalyser.jl LVg18g59 cshow
julia GlasgowPB3trimnalyser.jl LVg16g58 cshow
julia GlasgowPB3trimnalyser.jl LVg12g62 cshow

julia GlasgowPB3trimnalyser.jl LVg26g100 cshow
solve 14s
LVg26g100 & 63.85 MB &          & 4.682 MB & 0    & 0    & 40.6 ( 4.97 34.4 1.22 0 ) \\\hline


remarques, inj1 n'est pas utilise mais deg36=1 l'est
=#
function run_LV_solver()
    path = string(benchs,"/LV")
    cd()
    graphs = cd(readdir, path)
    n = length(graphs)
    if CONFIG.insid != ""
        ins = string("LV",CONFIG.insid)
        i = findlast(c -> c == 'g', CONFIG.insid)
        pattern = CONFIG.insid[1:i-1]
        target = CONFIG.insid[i:end]
        solve(ins,path,pattern,path,target,"lad",0,1_000_000_000)
        if CONFIG.veripb
            cd()
            cd(CONFIG.pbopath)
            runpboxide(ins)
        end
    else
    t=0
    stats = [stat(path*'/'*file).size for file in graphs]
    println("stats: ", stats)
    p = sortperm(stats)
    for i in eachindex(graphs), j in eachindex(graphs)
        pattern = graphs[p[i]]
        target = graphs[p[j]]
    # for target in graphs[1:end], pattern in graphs[1:end]
        if pattern != target
            ins = string("LV",pattern,target)
            solve(ins,path,pattern,path,target,"lad")
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
        end
    end end
end
function solve(ins,pathpat,pattern,pathtar,target,format,minsize=1_000,maxsize=50_000_000,remake=true,verbose=false)
    if remake || !isfile(string(proofs,"/",ins,".opb")) || !isfile(string(proofs,"/",ins,extention)) || 
            length(read(`tail -n 1 $proofs/$ins$extention`,String)) < 24 ||
            read(`tail -n 1 $proofs/$ins$extention`,String)[1:24] != "end pseudo-Boolean proof"
        print(ins,' ')
        cd(CONFIG.solveurpath)
        t = @elapsed begin
            # p = run(pipeline(`timeout $timeout ./$solver --prove $proofs/$ins --no-supplementals --no-clique-detection --format $format $pathpat/$pattern $pathtar/$target`, devnull),wait=false); wait(p)
            try
                if verbose
                    p = run(`timeout $tl ./glasgow_subgraph_solver --prove $proofs/$ins --no-clique-detection --format $format $pathpat/$pattern $pathtar/$target`)
                else
                    redirect_stdio(stdout = devnull,stderr = devnull) do
                    p = run(`timeout $tl ./glasgow_subgraph_solver --prove $proofs/$ins --no-clique-detection --format $format $pathpat/$pattern $pathtar/$target`)
                    end
                end
            catch e
            end
        end
        t+=0.01
        ok = false
        print(prettytime(t))
        if t>tl
            printstyled(" timeout "; color = :red)
        elseif read(`tail -n 2 $proofs/$ins$extention`,String)[1:14] == "conclusion SAT"
            printstyled(" sat     "; color = 166)
        elseif minsize > stat(string(proofs,"/",ins,extention)).size            
            printstyled(" toosmal "; color = :yellow)
        elseif maxsize < stat(string(proofs,"/",ins,extention)).size            
            printstyled(" toobig "; color = :red)
        else printstyled(" OK      "; color = :green)
            ok = true
            # g = ladtograph(pathpat,pattern)
            # draw(PNG(string(proofs,"/aimg/graphs/",ins,pattern[1:3],".png"), 16cm, 16cm), gplot(g))
            # g = ladtograph(pathtar,target)
            # draw(PNG(string(proofs,"/aimg/graphs/",ins,target[1:3],".png"), 16cm, 16cm), gplot(g))
        end
        println()
        if !ok
            run(`rm -f $proofs/$ins$extention`)
            run(`rm -f $proofs/$ins.opb`)
        end
    end
end
function runpboxide(file)
    t1 = t2 = 0
    prof = CONFIG.profiling ? "flamegraph" : "r"
    t1 = @elapsed begin
        try
            printstyled("veriPB check ", color=:blue)
                if CONFIG.trace
                    v1 = run(`timeout $tl cargo $prof -- $proofs/$file.opb $proofs/$file$extention `)
                else
                    redirect_stdio(stdout = devnull,stderr = devnull) do
                        v1 = run(`timeout $tl cargo $prof -- $proofs/$file.opb $proofs/$file$extention `)
                    end
                end
        catch e
            printstyled("fail ", color=:red)
            # print(e)
        end
    end
    printstyled(prettytime(t1),"  \n"; color = :cyan)
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
main()