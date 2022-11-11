export Brownian
export update_particles!

mutable struct Brownian
    particles::Vector{ActiveParticle}

    dt::Float64
end

function update_particles!(integrator::Brownian; box::Vector{Float64})
    @inbounds Threads.@threads for particle in integrator.particles
        if particle.run_time_remaining <= 0.0
            norm² = particle.preferred_director[1]^2 + particle.preferred_director[2]^2
            if norm² > 0.0
                norm = sqrt(norm²)
                sin, cos = sincos(particle.noise * randn())
                dir_x, dir_y = particle.preferred_director[1] / norm, particle.preferred_director[2] / norm

                particle.director[1] = cos * dir_x - sin * dir_y
                particle.director[2] = sin * dir_x + cos * dir_y
            else
                particle.director[2], particle.director[1] = sincos(2 * pi * rand())
            end
            particle.run_time_remaining = -particle.decision_time * log(rand())
        end

        particle.position[1] += integrator.dt * (particle.director[1] * particle.speed + particle.force[1] / particle.γ_trans)
        particle.position[2] += integrator.dt * (particle.director[2] * particle.speed + particle.force[2] / particle.γ_trans)

        if particle.position[1] < 0.0 || particle.position[1] > box[1]
            particle.position[1] = mod(particle.position[1], box[1])
        end
        if particle.position[2] < 0.0 || particle.position[2] > box[2]
            particle.position[2] = mod(particle.positionp[2], box[2])
        end

        fill!(particle.force, 0.0)
        fill!(particle.preferred_director, 0.0)
        particle.run_time_remaining -= integrator.dt
    end
end