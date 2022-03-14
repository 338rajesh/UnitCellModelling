
"""
    modelling() -> Dict

Returns a julia dictionary with the following key-value pairs

- `"mesh_data"` => `Dict{String, Any}` , where possible key-value pairs of mesh data are
    - `"all_node_tags"` => `Vector{Int}`
    - `"all_node_coordinates"` => `Matrix{Float64}`
    - `"matrix_element_connectivity"` => `Matrix{Int}`
    - `"inclusions_element_connectivity"` => `Matrix{Int}`
    - `"total_nodes"` => Int
    - `"num_2D_elements"` => Int
    - `"num_3D_elements"` => Int


"""
function create(
    unit_cell::AbstractUnitCell,
    options::Union{Tuple,Dict{Symbol,Any}} = ()
)::Dict{String,Any}
    model_data::Dict{String, Any} = Dict{String, Any}()
    unit_cell_dim = dimension(unit_cell)
    if isempty(options)  # FIXME update the way default options are set!
        options = Dict(
            :small_parameter => 1e-06,
            :geom_export_paths => (joinpath(homedir(), "unit_cell.step"),),
            :extrusion_direction => "XY->Z",
            :mesh_periodicity => true,
            :element_types => (:DEFAULT,),
            :num_ele_in_extr_dir => Int64[],
            :cum_heights => Float64[],
            :recombine_ele_in_extr_dir => true,
            :mesh_min_size_factor => 2.0,
            :mesh_max_size_factor => 1.0,
            :mesh_opt_algorithm => "Netgen",
            :show_mesh_stats => true,
            :show_rve => true,
        )
    end
    global ϵ = options[:small_parameter]
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
    ucp_geom_tags = geometric_model(unit_cell, options)  # unit cell phase tags
    matrix_pg_tag = gmsh.model.add_physical_group(unit_cell_dim, ucp_geom_tags[1:1])
    inclusions_pg_tags = gmsh.model.add_physical_group(unit_cell_dim, ucp_geom_tags[2:end])
    gmsh.model.set_physical_name(unit_cell_dim, matrix_pg_tag, "Matrix")
    gmsh.model.set_physical_name(unit_cell_dim, inclusions_pg_tags, "Inclusions")
    #
    # ---------------------
    #       Meshing
    # ---------------------
    generate_mesh(unit_cell, options)
    export_geometry(options)
    #
    # ---------------------
    #       EXPort
    # ---------------------
    model_data["mesh_data"] = get_mesh_data(unit_cell, ucp_geom_tags)  # Dict{String, Any}
    if options[:show_rve]
        gmsh.fltk.run()
    end
    #
    gmsh.finalize()
    #
    return model_data
end
