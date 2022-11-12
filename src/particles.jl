mutable struct ActiveParticle
    position::MVector{2, Float64}
    director::MVector{2, Float64}

    force::MVector{2, Float64}
    preferred_director::MVector{2, Float64}

    type::Symbol
    speed::Float64
    decision_time::Float64
    run_time_remaining::Float64
    noise::Float64

    # constants for now
    γ_trans::Float64
    R::Float64
end

function ActiveParticle(type::Symbol, position::Vector{Float64}, speed::Float64, decision_time::Float64; noise::Float64 = 0.0, γ_trans::Float64 = 1.0, R::Float64 = 0.5)
    run_time_remaining = -decision_time * log(rand())
    sin, cos = sincos(2 * pi * rand())
    return ActiveParticle(position, [cos, sin], zeros(2), zeros(2), type, speed, decision_time, run_time_remaining, noise, γ_trans, R)
end

function group_by_type(particles::Vector{ActiveParticle}, types::Vector{Symbol})
    subset = Vector{ActiveParticle}()

    for particle in particles
        if particle.type in types
            push!(subset, particle) 
        end
    end

    return subset
end