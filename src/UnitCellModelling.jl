module UnitCellModelling

    using gmsh
    # using 

    include("ucmBase.jl")
    include("gshapes.jl")
    include("geometric_modelling.jl")

    export model

end
