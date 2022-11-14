using Plots

export visualize!

function visualize!(system::System; save_as::String = "TEMP/TEMP.gif", fps::Int64 = 10, frame_nums::Vector{Int64} = Int64[])
    if !isdir(dirname(save_as))
        mkpath(dirname(save_as))
    end

    if isempty(frame_nums)
        frame_nums = [f for f = 1 : length(system.history)]
    end

    xlim, ylim = [0.0, system.box[1]], [0.0, system.box[2]]

    println("GENERATING GIF")
    animation = Animation()

    N_particles = length(system.particles)
    scene = plot(size = (600, 600), legend = false, axis = false, grid = false, aspect_ratio = 1, xlim = xlim, ylim = ylim)
    for (f, frame_num) = enumerate(frame_nums)
        if f == 1
            θ = LinRange(0.0, 2 * pi, 20)
            unit_circle = [cos.(θ), sin.(θ)]
            for particle in system.history[frame_num]
                xs, ys = particle.R .* unit_circle[1] .+ particle.position[1], particle.R .* unit_circle[2] .+ particle.position[2]
                plot!(xs, ys, seriestype = [:shape,], color = :red, linecolor = :red, fillalpha = 0.3)
            end
        else
            prev_frame_num = frame_nums[f - 1]
            for n = 1 : N_particles
                prev_particle, particle = system.history[prev_frame_num][n], system.history[frame_num][n]
                Δx, Δy = particle.position[1] - prev_particle.position[1], particle.position[2] - prev_particle.position[2]
                scene.series_list[n][:x] .+= Δx
                scene.series_list[n][:y] .+= Δy
            end
        end

        frame(animation)
        println("Frame ", frame_num, " done")
    end
    gif(animation, save_as, fps = fps)
end