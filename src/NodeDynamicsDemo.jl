module NodeDynamicsDemo

using Plots


struct EnergyNode
    energyContent::Vector{Float64}
    selfDischarge::Float64
    diffCoeffPrev::Float64
    diffCoeffNext::Float64
    influx::Vector{Float64}
end


struct NodeChain
    nodes::Vector{EnergyNode}
    fix_initial_boundary::Vector{Float64}
    fix_final_boundary::Vector{Float64}
    timesteps::Vector{Float64}
    function NodeChain(
        nodes::Vector{EnergyNode},
        fix_initial_boundary::Vector{Float64},
        fix_final_boundary::Vector{Float64},
        timesteps::Vector{Float64},
    )
        # Check that all time series have proper lengths.
        lenvec = length.([n.influx for n in nodes])
        push!(lenvec, length(fix_initial_boundary))
        push!(lenvec, length(fix_final_boundary))
        if any(length(timesteps) .!= lenvec)
            @error "The provided vector data doesn't match in length!"
        end
        new(
            nodes,
            fix_initial_boundary,
            fix_final_boundary,
            timesteps,
        )
    end
end


function get_node_energy(c::NodeChain, j::Integer, i::Integer, f::Symbol)
    n = get(c.nodes, j, nothing)
    if isnothing(n)
        return getfield(c, f)[i]
    else
        return n.energyContent[i]
    end
end


function solve_explicit_dynamics!(chain::NodeChain)
    # Loop over the timesteps
    for (i, dt) in enumerate(chain.timesteps)
        # Loop over the nodes
        for (j, n) in enumerate(chain.nodes)
            e_prev = last(n.energyContent)
            e_next = e_prev + dt * (
                n.diffCoeffPrev * (
                    get_node_energy(chain, j - 1, i, :fix_initial_boundary) - e_prev
                )
                + n.diffCoeffNext * (
                    get_node_energy(chain, j + 1, i, :fix_final_boundary) - e_prev
                )
                + n.influx[i]
                -
                n.selfDischarge * e_prev
            )
            push!(n.energyContent, e_next)
        end
    end
end


function plot_chain(chain::NodeChain; kwargs...)
    plot(
        cumsum(vcat([0], chain.timesteps)), # Add initial value time step at zero
        hcat(getfield.(chain.nodes, :energyContent)...);
        kwargs...
    )
end

export EnergyNode, NodeChain, solve_explicit_dynamics!, plot_chain

end # module NodeDynamicsDemo
