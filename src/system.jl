mutable struct System
    descriptor::String
    
    box::Vector{Float64}

    particles::Vector{ActiveParticle}
    neighbor_lists::Vector{LinkedCellList}
    interactions::Vector{AbstractInteraction}
    integrators::Vector{Brownian}

    history::Vector{Vector{ActiveParticle}}
end

function System()
    return System("", zeros(2), Vector{ActiveParticle}(), Vector{LinkedCellList}(), Vector{AbstractInteraction}(), Vector{Brownian}(), Vector{Vector{ActiveParticle}}())
end

function generate_lattice(unit_cell::Vector{Vector{Float64}}, lattice_vector::Vector{Vector{Float64}}, duplicate::Vector{Int64})
    L_x = lattice_vector[1][1] * duplicate[1]
    L_y = lattice_vector[2][2] * duplicate[2]

    positions = Vector{Vector{Float64}}()
    for i = 0 : duplicate[1] -1, j = 0 : duplicate[2] - 1
        for position in unit_cell
            x, y = position .+ lattice_vector[1] .* i .+ lattice_vector[2] .* j

            push!(positions, [mod(x, L_x), mod(y, L_y)])
        end
    end

    return positions, [L_x, L_y]
end

@inline function hr_min_sec(time::Float64)
    hours = trunc(Int64, time / 3600.0)
    minutes = trunc(Int64, mod(time, 3600.0) / 60.0)
    seconds = trunc(Int64, mod(time, 60.0))

    return string(hours < 10 ? "0" : "", hours, 
                  minutes < 10 ? ":0" : ":", minutes, 
                  seconds < 10 ? ":0" : ":", seconds)
end

function run_simulation!(system::System; num_steps::Int64 = 1, save_interval::Int64 = -1, message_interval::Float64 = 10.0)
    println("")
    println("---------------SIMULATION STARTED---------------")
    println("Number of particles: $(length(system.particles))")
    println("Box dimensions: $(round.(system.box; digits = 1))")
    println("Threads available: $(Threads.nthreads())")

    prev_step = 0
    time_elapsed = 0.0
    interval_start = time()
    @time for step = 1 : num_steps
        for interaction in system.interactions
            compute_interaction!(interaction; box = system.box)
        end

        if step % save_interval == 0 && save_interval > 0
            push!(system.history, deepcopy(system.particles))
        end

        for integrator in system.integrators
            update_particles!(integrator; box = system.box)
        end

        for neighbor_list in system.neighbor_lists
            if num_steps % neighbor_list.update_interval == 0
                update_neighbor_list!(neighbor_list)
            end
        end

        interval_time = time() - interval_start
        if interval_time > message_interval || step == num_steps
            time_elapsed += interval_time
            rate = (step - prev_step) / interval_time
            println(hr_min_sec(time_elapsed), " | ",
                    step, "/", num_steps, " (", round(step / num_steps * 100, digits = 1), "%) | ",
                    round(rate, digits = 1), " steps/s | ",
                    hr_min_sec((num_steps - step) / rate))
            prev_step = step
            interval_start = time()
        end
    end
    println("Average steps/s: ", round(num_steps / time_elapsed, digits = 1))
    println("----------------SIMULATION DONE-----------------")
    println("")
end

function save_simulation!(system::System; save_as::String = "TEMP.out")
    if !isdir(dirname(save_as))
        mkpath(dirname(save_as))
    end

    println("Saving simulation...")
    open(save_as, "w") do f
        serialize(f, system)
    end
    println("Simulation saved to ", save_as)
end

function load_simulation()
    system = begin
        open(file, "r") do f
            deserialize(f)
        end
    end

    return system
end

function generate_mathematica_data!()
    
end