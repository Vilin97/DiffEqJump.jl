using DiffEqJump, DiffEqBase

# SIR model
reactstoch = [
    [1 => 1, 2 => 1],
    [2 => 1],
]
netstoch = [
    [1 => -1, 2 => 1],
    [2 => -1, 3 => 1],
]
spec_to_dep_jumps = [[1],[1,2],convert(Array{Int64,1}, [])]
jump_to_dep_specs = [[1,2],[2,3]]
rates = [1e-4, 0.01]
majumps = MassActionJump(rates, reactstoch, netstoch)
prob = DiscreteProblem([999,1,0],(0.0,250.0), rates)
jump_prob_SIR = JumpProblem(prob, DirectCR(), majumps, save_positions=(false,false), vartojumps_map=spec_to_dep_jumps, jumptovars_map=jump_to_dep_specs)

integrator = init(jump_prob_SIR, SSAStepper())

using BenchmarkTools, Parameters
p = jump_prob_SIR.discrete_jump_aggregation
@unpack ma_jumps, rates, rng = p
@unpack t = integrator

#
# using Plots; plotly(); plot(sol)
#
# nums = Int[]
# @time for i in 1:1000
#   sol = solve(jump_prob,FunctionMap())
#   push!(nums,sol[end][1])
# end
# mean(nums)
#
# using ProfileView
# @profile for i in 1:1000; solve(jump_prob,FunctionMap()); end
# Profile.clear()
# @profile for i in 1:1000; solve(jump_prob,FunctionMap()); end
# ProfileView.view()
