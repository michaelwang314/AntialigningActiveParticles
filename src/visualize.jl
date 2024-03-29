using Plots

export visualize!

function degrees(director::V) where V <: AbstractVector
    θ = atan(director[2], director[1])
    return (θ + (θ < 0.0 ? 2 * pi : 0.0)) * 180 / pi
end

function visualize!(system::System; save_as::String = "TEMP/TEMP.gif", fps::Int64 = 10, frame_nums::Vector{Int64} = Int64[], colors::Dict{Symbol, Symbol} = Dict{Symbol,Symbol}(:notype => :black))
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
            unit_circle = [cos.(LinRange(0.0, 2 * pi, 20)), sin.(LinRange(0.0, 2 * pi, 20))]
            for particle in system.history[frame_num]
                xs, ys = particle.R .* unit_circle[1] .+ particle.position[1], particle.R .* unit_circle[2] .+ particle.position[2]
                
                if colors[particle.type] == :wheel
                    θ = degrees(particle.director)
                    color = RGB(HSV(θ, 1.0, 1.0))
                else
                    color = colors[particle.type]
                end
                
                plot!(xs, ys, seriestype = [:shape,], fillcolor = color, linecolor = color, fillalpha = 0.3)
            end
        else
            prev_frame_num = frame_nums[f - 1]
            for n = 1 : N_particles
                prev_particle, particle = system.history[prev_frame_num][n], system.history[frame_num][n]
                Δx, Δy = particle.position[1] - prev_particle.position[1], particle.position[2] - prev_particle.position[2]
                scene.series_list[n][:x] .+= Δx
                scene.series_list[n][:y] .+= Δy

                if colors[particle.type] == :wheel
                    θ = degrees(particle.director)
                    color = RGB(HSV(θ, 1.0, 1.0))

                    scene.series_list[n][:fillcolor] = color
                    scene.series_list[n][:linecolor] = color
                end
            end
        end

        frame(animation)
        println("Frame ", frame_num, " done")
    end
    gif(animation, save_as, fps = fps)
end