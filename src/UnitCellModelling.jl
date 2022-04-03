module UnitCellModelling

using gmsh
using Parameters
using StaticArrays
# using 

include("ucmBase.jl")
include("gshapes.jl")
include("geometric_modelling.jl")
include("meshing.jl")

"""
    make_unit_cell_model() -> Dict


It performs mainly three tasks
 1. Geometry creation
 2. Mesh generation
 3. Exporting geometry and mesh data

Returns a julia dictionary with the following key-value pairs

- `"mesh_data"` => `Dict{String, Any}` , where possible key-value pairs of dictionary are
    - `"all_node_tags"` => `Vector{Int}`
    - `"all_node_coordinates"` => `Matrix{Float64}`
    - `"matrix_element_connectivity"` => `Matrix{Int}`
    - `"inclusions_element_connectivity"` => `Matrix{Int}`
    - `"total_nodes"` => Int
    - `"num_2D_elements"` => Int
    - `"num_3D_elements"` => Int


"""
function make_unit_cell_model(
    unit_cell::AbstractUnitCell;
    small_parameter::Float64 = 1e-06,
    geom_export_paths::Tuple{Vararg{String}} = (),  # joinpath(homedir(), "unit_cell.step"),
    extr_dir::String = "XY->Z",
    extr_dir_num_ele::Vector{Int64} = Int64[],
    extr_dir_cum_heights::Vector{Float64} = Float64[],
    extr_dir_recombine_ele::Bool = true,
    mesh_periodicity::Bool = true,
    element_types::Tuple{Vararg{Symbol}} = (:DEFAULT,),
    min_ele_size_factor::Float64 = 1.0,   # FIXME make default to take factor, instead of absolute values
    max_ele_size_factor::Float64 = 2.0,
    mesh_opt_algorithm::String = "Netgen",
    show_mesh_stats::Bool = true,
    show_rve::Bool = true,
    verbose::Int=1,
)::Dict{String,Any}
    unit_cell_dim::Int = dimension(unit_cell)
    small_uc_side_length = minimum(side_lengths(unit_cell))
    global ϵ = small_parameter

    #
    gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", 0)
    gmsh.option.set_number("General.Terminal", 0)
    gmsh.option.set_number("General.NoPopup", 0)
    gmsh.option.set_number("Geometry.OCCBoundsUseStl", 1)
    gmsh.model.add("unit_cell")
    #
    # ---------------------
    #   Geometry creation
    # ---------------------    
    ucp_geom_tags = geometric_model(
        unit_cell,
        extr_dir,
        extr_dir_num_ele,
        extr_dir_cum_heights,
        extr_dir_recombine_ele
    )  # unit cell phase tags
    matrix_pg_tag = gmsh.model.add_physical_group(unit_cell_dim, ucp_geom_tags[1:1])
    inclusions_pg_tags = gmsh.model.add_physical_group(unit_cell_dim, ucp_geom_tags[2:end])
    gmsh.model.set_physical_name(unit_cell_dim, matrix_pg_tag, "Matrix")
    gmsh.model.set_physical_name(unit_cell_dim, inclusions_pg_tags, "Inclusions")
    #
    # ---------------------
    #       Meshing
    # ---------------------
    generate_mesh(
        unit_cell,
        mesh_periodicity,
        element_types,
        extr_dir_num_ele,
        min_ele_size_factor*small_uc_side_length,
        max_ele_size_factor*small_uc_side_length,
        mesh_opt_algorithm,
    )
    if show_mesh_stats && (verbose > 0)
        write_mesh_statistics()
    end
    #
    # ---------------------
    #       EXPort
    # ---------------------

    # exporting UC model to the specified fomats
    if !isempty(geom_export_paths)
        for a_export_path in geom_export_paths
            gmsh.write(a_export_path)
        end
    end

    # 
    model_data::Dict{String,Any} = Dict{String,Any}()
    model_data["mesh_data"] = get_mesh_data(unit_cell_dim, ucp_geom_tags)  # Dict{String, Any}
    #
    if show_rve
        gmsh.fltk.run()
    end
    #
    gmsh.finalize()
    #
    return model_data
end  # of function

export make_unit_cell_model

end  # of UnitCellModelling module
