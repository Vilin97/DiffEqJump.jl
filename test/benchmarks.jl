using DiffEqBase, DiffEqJump
using BenchmarkTools

alg = RDirect()
counter_coeffs = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0]

function time_Nsims(jprob, Nsims)
    solve(jump_prob, SSAStepper())
    res = 0
    for i in 1:Nsims
        res += @elapsed solve(jump_prob, SSAStepper())
    end
    res
end

# bimolerx
tf           = .01
u0           = [200, 100, 150]
reactstoch = [
    [1 => 2],
    [2 => 1],
    [1 => 1, 2 => 1],
    [3 => 1],
    [3 => 3]
]
netstoch = [
    [1 => -2, 2 => 1],
    [1 => 2, 2 => -1],
    [1 => -1, 2 => -1, 3 => 1],
    [1 => 1, 2 => 1, 3 => -1],
    [1 => 3, 3 => -3]
]
rates = [1., 2., .5, .75, .25]
spec_to_dep_jumps = [[1,3],[2,3],[4,5]]
jump_to_dep_specs = [[1,2],[1,2],[1,2,3],[1,2,3],[1,3]]
majumps = MassActionJump(rates, reactstoch, netstoch)
prob = DiscreteProblem(u0, (0.0, tf), rates)

println("Timing bimolerx")
for counter_coeff in counter_coeffs
    counter_threshold = trunc.(Int, counter_coeff*length(rates))
    jump_prob = JumpProblem(prob, alg, majumps, vartojumps_map=spec_to_dep_jumps, jumptovars_map=jump_to_dep_specs, counter_threshold = counter_threshold)
    t = time_Nsims(jump_prob, 1000)
    @show counter_coeff, t
end

# extinction
reactstoch = [
    [1 => 1]
]
netstoch = [
    [1 => -1]
]
rates = [1.]
spec_to_dep_jumps = [[1]]
jump_to_dep_specs = [[1]]
dg = [[1]]
majumps = MassActionJump(rates, reactstoch, netstoch)
u0 = [10]
prob = DiscreteProblem(u0,(0.,100.),rates)

println("Timing extinction")
ts = zeros(length(counter_coeffs))
for (i,counter_coeff) in enumerate(counter_coeffs)
    counter_threshold = trunc.(Int, counter_coeff*length(rates))
    jump_prob = JumpProblem(prob, alg, majumps, vartojumps_map=spec_to_dep_jumps, jumptovars_map=jump_to_dep_specs, counter_threshold = counter_threshold)
    t = time_Nsims(jump_prob, 1000)
    ts[i] = t
    @show counter_coeff, t
end

# BCR
using ReactionNetworkImporters, DiffEqBiological, OrdinaryDiffEq
# networkname = "tester"
# tf = 0.001
# # bngfname = "PATH/TO/FILE"
# bngfname = joinpath(dirname(pathof(ReactionNetworkImporters)),"..","data","bcr","bcr.net")
# prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), bngfname);
# rn = prnbng.rn; u₀ = prnbng.u₀; p = prnbng.p; shortsymstosyms = prnbng.symstonames;
# u0 = round.(Int, u₀)
# addjumps!(rn,build_regular_jumps=false, minimal_jumps=true)
# prob = DiscreteProblem(rn, u0, (0.,tf), p)
#
# ts = zeros(length(counter_coeffs))
# for (i,counter_coeff) in enumerate(counter_coeffs)
#     counter_threshold = trunc.(Int,counter_coeff*24388)
#     jump_prob = JumpProblem(prob, RDirect(), rn, save_positions=(false,false), counter_threshold = counter_threshold)
#     t = time_Nsims(jump_prob, 10)
#     ts[i] = t
#     @show counter_coeff, t
# end
# println(ts)

#diffusion
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

prob = DiscreteProblem(u0, (0.0, tf), rnpar)

ts = zeros(length(counter_coeffs))
for (i,counter_coeff) in enumerate(counter_coeffs)
    counter_threshold = trunc.(Int,counter_coeff*2*N)
    jump_prob = JumpProblem(prob, RDirect(), rn, save_positions=(false,false), counter_threshold = counter_threshold)
    t = time_Nsims(jump_prob, 100)
    ts[i] = t
    @show counter_coeff, t
end
println(ts)

# gene expression
reactstoch = [
    [1 => 1],
    [2 => 1],
    [2 => 1],
    [3 => 1],
    [1 => 1, 3 => 1],
    [4 => 1]
]
netstoch = [
    [2 => 1],
    [3 => 1],
    [2 => -1],
    [3 => -1],
    [1 => -1, 3 => -1, 4 => 1],
    [1 => 1, 3 => 1, 4 => -1]
]
spec_to_dep_jumps = [[1,5],[2,3],[4,5],[6]]
jump_to_dep_specs = [[2],[3],[2],[3],[1,3,4],[1,3,4]]
rates = [.5, (20*log(2.)/120.), (log(2.)/120.), (log(2.)/600.), .025, 1.]
majumps = MassActionJump(rates, reactstoch, netstoch)
u0           = [1,0,0,0]
tf           = 1000.0
# TESTING:
prob = DiscreteProblem(u0, (0.0, tf), rates)
ts = zeros(length(counter_coeffs))
for (i,counter_coeff) in enumerate(counter_coeffs)
    counter_threshold = trunc.(Int,counter_coeff*6)
    jump_prob = JumpProblem(prob, RDirect(), majumps, save_positions=(false,false), counter_threshold = counter_threshold)
    t = time_Nsims(jump_prob, 100000)
    ts[i] = t
    @show counter_coeff, t
end
jump_prob = JumpProblem(prob, RDirect(), majumps, save_positions=(false,false), counter_threshold = 1)
