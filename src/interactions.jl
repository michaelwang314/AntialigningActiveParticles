mutable struct LennardJones
    ϵ::Float64
    σ::Float64
    cutoff::Float64

    particles::Vector{ActiveParticle}
    neighbor_list::LinkedCellList
end

mutable struct Alignment
    align::Float64
    cutoff::Float64

    particles::Vector{ActiveParticle}
    neighbor_list::LinkedCellList
end

mutable struct CircularConfinement
    particles::Vector{ActiveParticle}

    ϵ::Float64
    σ::Float64
    radius::Float64
    center::Vector{Float64}
end

function compute_interaction!(interaction::LennardJones; box::Vector{Float64})

end

function compute_interaction!(interaction::Alignment; box::Vector{Float64})

end

function compute_interaction!(interaction::CircularConfinement; args...)
    
end