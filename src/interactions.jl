abstract type AbstractInteraction end

mutable struct LennardJones <: AbstractInteraction
    ϵ::Float64
    σ::Float64
    cutoff::Float64

    particles::Vector{ActiveParticle}
    neighbor_list::LinkedCellList
end

mutable struct Alignment <: AbstractInteraction
    angle::Float64
    cutoff::Float64

    particles::Vector{ActiveParticle}
    neighbor_list::LinkedCellList
end

mutable struct CircularConfinement <: AbstractInteraction
    particles::Vector{ActiveParticle}

    ϵ::Float64
    σ::Float64
    radius::Float64
    center::Vector{Float64}
    side::Symbol
end

mutable struct ChannelConfinement <: AbstractInteraction
    particles::Vector{ActiveParticle}

    ϵ::Float64
    σ::Float64
    width::Float64
    center_y::Float64
end

mutable struct Gravity <: AbstractInteraction
    particles::Vector{ActiveParticle}

    mg::Float64
end

function compute_interaction!(lennardjones::LennardJones; box::Vector{Float64})
    @inbounds Threads.@threads for particle in lennardjones.particles
        x, y = particle.position
        i = trunc(Int64, x / lennardjones.neighbor_list.cell_spacing_x)
        j = trunc(Int64, y / lennardjones.neighbor_list.cell_spacing_y)
        for di = -1 : 1, dj = -1 : 1
            idi = mod(i + di, lennardjones.neighbor_list.num_cells_x) + 1
            jdj = mod(j + dj, lennardjones.neighbor_list.num_cells_y) + 1
            id = lennardjones.neighbor_list.start_id[idi, jdj]
            while id > 0
                neighbor = lennardjones.neighbor_list.particles[id]
                Δ = x - neighbor.position[1]
                Δx = Δ - (i + di + 1 != idi ? sign(Δ) * box[1] : 0.0)
                Δ = y - neighbor.position[2]
                Δy = Δ - (j + dj + 1 != jdj ? sign(Δ) * box[2] : 0.0)
                Δr² = Δx^2 + Δy^2

                if 0.0 < Δr² < lennardjones.cutoff^2
                    val = (lennardjones.σ^2 / Δr²)^3
                    coef = lennardjones.ϵ * (48.0 * val - 24.0) * val / Δr²

                    particle.force[1] += coef * Δx
                    particle.force[2] += coef * Δy
                end
                id = lennardjones.neighbor_list.next_id[id]
            end
        end
    end
end

function compute_interaction!(alignment::Alignment; box::Vector{Float64})
    @inbounds Threads.@threads for particle in alignment.particles
        if particle.run_time_remaining <= 0.0
            x, y = particle.position
            i = trunc(Int64, x / alignment.neighbor_list.cell_spacing_x)
            j = trunc(Int64, y / alignment.neighbor_list.cell_spacing_y)
            for di = -1 : 1, dj = -1 : 1
                idi = mod(i + di, alignment.neighbor_list.num_cells_x) + 1
                jdj = mod(j + dj, alignment.neighbor_list.num_cells_y) + 1
                id = alignment.neighbor_list.start_id[idi, jdj]
                while id > 0
                    neighbor = alignment.neighbor_list.particles[id]
                    Δ = x - neighbor.position[1]
                    Δx = Δ - (i + di + 1 != idi ? sign(Δ) * box[1] : 0.0)
                    Δ = y - neighbor.position[2]
                    Δy = Δ - (j + dj + 1 != jdj ? sign(Δ) * box[2] : 0.0)

                    if 0.0 < Δx^2 + Δy^2 < alignment.cutoff^2
                        ndirx, ndiry = neighbor.director
                        pdirx, pdiry = particle.director
                        sin, cos = sincos(sign((ndirx * pdirx + ndiry * pdiry) * (Δx * pdiry - pdirx * Δy)) * alignment.angle)
                        particle.preferred_director[1] += cos * ndirx - sin * ndiry
                        particle.preferred_director[2] += sin * ndirx + cos * ndiry
                    end
                    id = alignment.neighbor_list.next_id[id]
                end
            end
        end
    end
end

function compute_interaction!(circularconfinement::CircularConfinement; args...)
    side = circularconfinement.side == :in ? 1.0 : -1.0
    @inbounds Threads.@threads for particle in circularconfinement.particles
        rx = particle.position[1] - circularconfinement.center[1]
        ry = particle.position[2] - circularconfinement.center[2]
        r² = rx^2 + ry^2

        if side * ((circularconfinement.radius - side * circularconfinement.σ * 2^(1.0 / 6.0))^2 - r²) < 0
            r = sqrt(r²)
            Δr² = (circularconfinement.radius - r)^2
            val = (circularconfinement.σ^2 / Δr²)^3
            coef = circularconfinement.ϵ * (48.0 * val - 24.0) * val / sqrt(Δr²) * (-side / r)

            particle.force[1] += coef * rx
            particle.force[2] += coef * ry
        end
    end
end

function compute_interaction!(channelconfinement::ChannelConfinement; args...)
    @inbounds Threads.@threads for particle in channelconfinement.particles
        ry = particle.position[2] - channelconfinement.center_y
        
        if (channelconfinement.width / 2 - channelconfinement.σ * 2^(1.0 / 6.0))^2 < ry^2
            Δr² = (channelconfinement.width / 2 - abs(ry))^2
            val = (channelconfinement.σ^2 / Δr²)^3
            coef = channelconfinement.ϵ * (48.0 * val - 24.0) * val / sqrt(Δr²)

            particle.force[2] += coef * sign(-ry)
        end
    end
end

function compute_interaction!(gravity::Gravity; args...)
    @inbounds Threads.@threads for particle in gravity.particles
        particle.force[2] -= gravity.mg
    end
end