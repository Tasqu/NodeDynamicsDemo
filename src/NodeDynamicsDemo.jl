module NodeDynamicsDemo

using Plots
using LinearAlgebra


struct EnergyNode
    energyContent::Vector{Float64}
    selfDischarge::Float64
    diffCoeffPrev::Float64
    diffCoeffNext::Float64
    influx::Vector{Float64}
    function EnergyNode(
        energyContent::Vector{Float64},
        selfDischarge::Float64,
        diffCoeffPrev::Float64,
        diffCoeffNext::Float64,
        influx::Vector{Float64}
    )
        # Check initial energyContent and sizehint
        if length(energyContent) != 1
            @error "`energyContent` should only have one element representing the initial value!"
        else
            sizehint!(energyContent, length(influx))
            new(
                energyContent,
                selfDischarge,
                diffCoeffPrev,
                diffCoeffNext,
                influx
            )
        end
    end
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


function solve_implicit_dynamics!(chain::NodeChain)
    # Form the inverse dynamics matrices for unique timestep lenghts
    invmats = Dict(
        dt => inv(
            Tridiagonal(
                [-n.diffCoeffPrev for n in chain.nodes[2:end]],
                [
                    1 / dt +
                    n.selfDischarge +
                    n.diffCoeffPrev +
                    n.diffCoeffNext
                    for n in chain.nodes
                ],
                [-n.diffCoeffNext for n in chain.nodes[1:end-1]]
            )
        )
        for dt in unique(chain.timesteps)
    )
    # Loop over the timesteps and solve the states.
    for (i, dt) in enumerate(chain.timesteps)
        # Calculate the initial state vector.
        E_prev = [last(n.energyContent) / dt + n.influx[i] for n in chain.nodes]
        E_prev[1] += first(chain.nodes).diffCoeffPrev * chain.fix_initial_boundary[i]
        E_prev[end] += last(chain.nodes).diffCoeffNext * chain.fix_final_boundary[i]
        # Solve the dynamics and save the results
        E_next = invmats[dt] * E_prev
        push!.(getfield.(chain.nodes, :energyContent), E_next)
    end
end


function solve_analytical_dynamics!(chain::NodeChain)
    # Form the dynamics matrix, its inverse, and exponentials.
    dynmat = Tridiagonal(
        [n.diffCoeffPrev for n in chain.nodes[2:end]],
        [-n.selfDischarge - n.diffCoeffPrev - n.diffCoeffNext for n in chain.nodes],
        [n.diffCoeffNext for n in chain.nodes[1:end-1]],
    )
    invmat = inv(dynmat)
    expmats = Dict(
        dt => exp(dynmat .* dt)
        for dt in chain.timesteps
    )
    # Loop over the timesteps and solve the states
    for (i, dt) in enumerate(chain.timesteps)
        # Fetch current states and influxes
        E_prev = [last(n.energyContent) for n in chain.nodes]
        P_prev = [n.influx[i] for n in chain.nodes]
        P_prev[1] += first(chain.nodes).diffCoeffPrev * chain.fix_initial_boundary[i]
        P_prev[end] += last(chain.nodes).diffCoeffNext * chain.fix_final_boundary[i]
        # Calculate and save the the next values
        E_next = expmats[dt] * E_prev - (I - expmats[dt]) * invmat * P_prev
        push!.(getfield.(chain.nodes, :energyContent), E_next)
    end
end


function plot_chain(chain::NodeChain; kwargs...)
    plot(
        cumsum(vcat([0], chain.timesteps)), # Add initial value time step at zero
        hcat(getfield.(chain.nodes, :energyContent)...);
        kwargs...
    )
end

export EnergyNode, NodeChain, solve_explicit_dynamics!, solve_implicit_dynamics!,
    solve_analytical_dynamics!, plot_chain

end # module NodeDynamicsDemo
