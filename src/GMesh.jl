"""
# GMesh

GMesh module is intended to do geometric modelling and meshing of Unit Cell using `gmsh`, an open source code.

It is designed to return the following information
1. Node Tags/Labels and corresponding coordinates
2. Element Tags/Labels and corresponding nodal_connectivity

"""
module GMesh
    using Reexport
    using gmsh, StaticArrays
    using ..UCHbase

    include("gshapes.jl")
    include("geometric_modelling.jl")
    include("meshing.jl")
    include("EXPort.jl")
    include("modelling.jl")
end
