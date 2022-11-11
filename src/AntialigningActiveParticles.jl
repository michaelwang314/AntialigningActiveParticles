module AntialigningActiveParticles
    using Serialization
    using StaticArrays
    using Plots

    export ActiveParticle, group_by_type
    export LinkedCellList
    export AbstractInteraction, LennardJones, Alignment, CircularConfinement
    export Brownian
    export System, generate_lattice, run_simulation!, save_simulation!, load_simulation!, generate_mathematica_data!


    """
    Flush output so that jobs can be monitored on a cluster.
    """
    @inline println(args...) = println(stdout, args...)
    @inline function println(io::IO, args...)
        Base.println(io, args...)
        flush(io)
    end

    include("particles.jl")
    include("neighborlist.jl")
    include("interactions.jl")
    include("integrators.jl")
    include("system.jl")
    
    include("visualize.jl")
end
