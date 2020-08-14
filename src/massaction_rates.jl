###############################################################################
# Stochiometry for a given reaction is a vector of pairs mapping species id to
# stochiometric coefficient.
###############################################################################

@inline @fastmath function evalrxrate(speciesvec::AbstractVector{T}, rxidx::S,
                              majump::MassActionJump{U,V,W})::R where {T,S,R,U <: AbstractVector{R},V,W}
    val = one(T)
    @inbounds for specstoch in majump.reactant_stoch[rxidx]
        specpop = speciesvec[specstoch[1]]
        val    *= specpop
        @inbounds for k = 2:specstoch[2]
            specpop -= one(specpop)
            val     *= specpop
        end
    end

    @inbounds return val * majump.scaled_rates[rxidx]
end

@inline @fastmath function executerx!(speciesvec::AbstractVector{T}, rxidx::S,
                                      majump::MassActionJump{U,V,W}) where {T,S,U,V,W}
    @inbounds net_stoch = majump.net_stoch[rxidx]
    @inbounds for specstoch in net_stoch
        speciesvec[specstoch[1]] += specstoch[2]
    end
    nothing
end

@inline @fastmath function executerx(speciesvec::SVector{T}, rxidx::S,
                                      majump::MassActionJump{U,V,W}) where {T,S,U,V,W}
    @inbounds net_stoch = majump.net_stoch[rxidx]
    @inbounds for specstoch in net_stoch
        speciesvec = setindex(speciesvec,speciesvec[specstoch[1]]+specstoch[2],specstoch[1])
    end
    speciesvec

    #=
    map(net_stoch) do stoch
        @inbounds speciesvec[stoch[1]] + stoch[2]
    end
    =#
end

function scalerates!(unscaled_rates::AbstractVector{U}, stochmat::AbstractVector{V}) where {U,S,T,W <: Pair{S,T}, V <: AbstractVector{W}}
    @inbounds for i in eachindex(unscaled_rates)
        coef = one(T)
        @inbounds for specstoch in stochmat[i]
            coef *= factorial(specstoch[2])
        end
        unscaled_rates[i] /= coef
    end
    nothing
end

function scalerate(unscaled_rate::U, stochmat::AbstractVector{Pair{S,T}}) where {U <: Number, S, T}
    coef = one(T)
    @inbounds for specstoch in stochmat
        coef *= factorial(specstoch[2])
    end
    unscaled_rate /= coef
end


###############################################################################
# dependency graph when MassActionJump uses pairs to represent (species,stoich)
###############################################################################

# map from species to Set of reactions depending on that species
function spec_to_dep_rxs_map(numspec, ma_jumps::MassActionJump)

    numrxs = get_num_majumps(ma_jumps)

    # map from a species to reactions that depend on it
    spec_to_dep_rxs = [Set{Int}() for n = 1:numspec]
    for rx in 1:numrxs
        for (spec,stoch) in ma_jumps.reactant_stoch[rx]
            push!(spec_to_dep_rxs[spec], rx)
        end
    end

    spec_to_dep_rxs
end

"make a map from reactions to dependent species"
function rxs_to_dep_spec_map(majumps)
    [[s for (s, c) in majumps.net_stoch[i]] for i in 1:get_num_majumps(majumps)]
end

# dependency graph is a map from a reaction to a vector of reactions
# that should depend on species it changes
function make_dependency_graph(numspec, ma_jumps::MassActionJump)

    numrxs          = get_num_majumps(ma_jumps)
    spec_to_dep_rxs = spec_to_dep_rxs_map(numspec, ma_jumps)

    # create map from rx to reactions depending on it
    dep_sets = [SortedSet{Int}() for n = 1:numrxs]
    for rx in 1:numrxs

        # rx changes spec, hence rxs depending on spec depend on rx
        for (spec,stoch) in ma_jumps.net_stoch[rx]
            for dependent_rx in spec_to_dep_rxs[spec]
                push!(dep_sets[rx], dependent_rx)
            end
        end
    end

    # convert to Vectors of Vectors
    dep_graph = Vector{Vector{Int}}(undef, numrxs)
    for rx = 1:numrxs
        dep_graph[rx] = [dep for dep in dep_sets[rx]]
    end

    dep_graph
end


# update dependency graph to make sure jumps depend on themselves
function add_self_dependencies!(dg)
    for (i,jump_deps) in enumerate(dg)
        if !any(y->isequal(y,i), jump_deps)
            push!(jump_deps, i)
            sort!(jump_deps)
        end
    end
end
