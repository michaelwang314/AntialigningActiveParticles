export LinkedCellList
export update_neighbor_list!

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

    return CellList(particles, start_id, next_id, num_cells_x, num_cells_y, cell_spacing_x, cell_spacing_y, update_interval)
end

function update_cell_list!(cell_list::CellList)
    fill!(cell_list.start_id, -1)
    fill!(cell_list.next_id, -1)

    @inbounds for (n, particle) in enumerate(cell_list.particles)
        i = trunc(Int64, particle.position[1] / cell_list.cell_spacing_x) + 1
        j = trunc(Int64, particle.position[2] / cell_list.cell_spacing_y) + 1

        if cell_list.start_id[i, j] > 0
            cell_list.next_id[n] = cell_list.start_id[i, j]
        end
        cell_list.start_id[i, j] = n
    end
end