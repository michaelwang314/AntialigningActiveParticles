struct LinkedCellList
    particles::Vector{ActiveParticle}
    start_id::Array{Int64, 2}
    next_id::Vector{Int64}

    num_cells_x::Int64
    num_cells_y::Int64
    cell_spacing_x::Float64
    cell_spacing_y::Float64

    update_interval::Int64
end

function LinkedCellList(particles::Vector{ActiveParticle}, box::Vector{Float64}, spacing::Float64; update_interval::Int64 = 1)
    num_cells_x = trunc(Int64, box[1] / spacing)
    num_cells_y = trunc(Int64, box[2] / spacing)
    cell_spacing_x = box[1] / num_cells_x
    cell_spacing_y = box[2] / num_cells_y

    start_id = -ones(Int64, num_cells_x, num_cells_y)
    next_id = -ones(Int64, length(particles))

    for (n, particle) in enumerate(particles)
        i = trunc(Int64, particle.position[1] / cell_spacing_x) + 1
        j = trunc(Int64, particle.position[2] / cell_spacing_y) + 1

        if start_id[i, j] > 0
            next_id[n] = start_id[i, j]
        end
        start_id[i, j] = n
    end

    return LinkedCellList(particles, start_id, next_id, num_cells_x, num_cells_y, cell_spacing_x, cell_spacing_y, update_interval)
end

function update_neighbor_list!(neighbor_list::LinkedCellList)
    fill!(neighbor_list.start_id, -1)
    fill!(neighbor_list.next_id, -1)

    @inbounds for (n, particle) in enumerate(neighbor_list.particles)
        i = trunc(Int64, particle.position[1] / neighbor_list.cell_spacing_x) + 1
        j = trunc(Int64, particle.position[2] / neighbor_list.cell_spacing_y) + 1

        if neighbor_list.start_id[i, j] > 0
            neighbor_list.next_id[n] = neighbor_list.start_id[i, j]
        end
        neighbor_list.start_id[i, j] = n
    end
end