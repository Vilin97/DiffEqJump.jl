using DiffEqBase, DiffEqBiological, DiffEqJump, Plots, Statistics, DataFrames
gr()
fmt = :png

N = 256
h = 1 / N
rn = @empty_reaction_network
function getDiffNetwork!(rn,N)
    for i = 1:N
        addspecies!(rn, Symbol(:u, i))
    end
    addparam!(rn, :β)
    for i = 1:N
        (i < N) && addreaction!(rn, :β, (Symbol(:u,i)=>1,), (Symbol(:u,i+1)=>1,))
        (i > 1) && addreaction!(rn, :β, (Symbol(:u,i)=>1,), (Symbol(:u,i-1)=>1,))
    end
    rn
end
getDiffNetwork!(rn,N)
addjumps!(rn, build_regular_jumps=false, minimal_jumps=true)
rnpar = [1/(h*h)]
u0 = 10*ones(Int64, N)
tf = .01

algs = DiffEqJump.JUMP_AGGREGATORS
shortlabels = [string(leg)[1:end-2] for leg in algs]
prob    = prob = DiscreteProblem(u0, (0.0, tf), rnpar)
# # ploth   = plot(reuse=false)
# for (i,method) in enumerate(algs)
#     println("Plotting method: ", method)
#     jump_prob = JumpProblem(prob, method, rn, save_positions=(false,false))
#     sol = solve(jump_prob, SSAStepper(), saveat=tf/1000.)
#     # plot!(ploth,sol.t,sol[Int(N//2),:],label=shortlabels[i], format=fmt)
# end
# # plot!(ploth, title="Population at middle lattice site", xlabel="time",format=fmt)

function run_benchmark!(t, jump_prob, stepper)
    sol = solve(jump_prob, stepper)
    @inbounds for i in 1:length(t)
        if i%10 == 1
            println("simulation $i starts")
        end
        t[i] = @elapsed (solve(jump_prob, stepper))
        if i%10 == 1
            println("simulation $i ended")
        end
    end
end

nsims = 50
benchmarks = Vector{Vector{Float64}}()
for method in algs
    println("Benchmarking method: ", method)
    jump_prob = JumpProblem(prob, method, rn, save_positions=(false,false))
    stepper = SSAStepper()
    t = Vector{Float64}(undef,nsims)
    run_benchmark!(t, jump_prob, stepper)
    push!(benchmarks, t)
end

medtimes = Vector{Float64}(undef,length(algs))
stdtimes = Vector{Float64}(undef,length(algs))
avgtimes = Vector{Float64}(undef,length(algs))
for i in 1:length(algs)
    medtimes[i] = median(benchmarks[i])
    avgtimes[i] = mean(benchmarks[i])
    stdtimes[i] = std(benchmarks[i])
end

println("setting up dataframe")
df = DataFrame(names=shortlabels,medtimes=medtimes,relmedtimes=(medtimes/medtimes[1]),
                avgtimes=avgtimes, std=stdtimes, cv=stdtimes./avgtimes)

                sa = [string(round(mt,digits=4),"s") for mt in df.medtimes]
bar(df.names,df.relmedtimes,legend=:false, fmt=fmt)
scatter!(df.names, .05 .+ df.relmedtimes, markeralpha=0, series_annotations=sa, fmt=fmt)
ylabel!("median relative to $(shortlabels[1])")
title!("256 Site 1D Diffusion CTRW")
