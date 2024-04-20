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
function make_udc_model(
    unit_cell::AbstractUnitCell;
    small_parameter::Float64=1e-06,
    geom_export_paths::Tuple{Vararg{String}}=(),  # joinpath(homedir(), "unit_cell.step"),
    extr_dir::String="XY->Z",
    extr_dir_num_ele::Vector{Int64}=Int64[],
    extr_dir_cum_heights::Vector{Float64}=Float64[],
    extr_dir_recombine_ele::Bool=true,
    mesh_periodicity::Bool=true,
    element_types::Tuple{Vararg{Symbol}}=(:DEFAULT,),
    min_ele_size_factor::Float64=1.0,   # FIXME make default to take factor, instead of absolute values
    max_ele_size_factor::Float64=2.0,
    mesh_opt_algorithm::String="Netgen",
    node_renum_algorithm::String="RCM",
    show_mesh_stats::Bool=true,
    show_rve::Bool=true,
    verbose::Int=1
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
        min_ele_size_factor * small_uc_side_length,
        max_ele_size_factor * small_uc_side_length,
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
    model_data["mesh_data"] = get_mesh_data(unit_cell_dim, ucp_geom_tags, node_renum_algorithm, verbose)
    model_data["side_lengths"] = side_lengths(unit_cell)
    #
    if show_rve
        gmsh.fltk.run()
    end
    #
    gmsh.finalize()
    #
    return model_data
end  # of function


function make_uc_geom_model(
    unit_cell::AbstractUnitCell;
    geom_export_paths::Tuple{Vararg{String}}=(),  # joinpath(homedir(), "unit_cell.step"),
    extr_dir::String="XY->Z",
    show_rve::Bool=true,
    verbose::Int=1
)
    global ϵ = small_parameter
    #
    gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", 0)
    gmsh.option.set_number("General.Terminal", 0)
    gmsh.option.set_number("General.NoPopup", 0)
    gmsh.option.set_number("Geometry.OCCBoundsUseStl", 1)
    gmsh.model.add("unit_cell")
    #
    # creating geometry
    verbose > 0 ? println("> creating geometric model") : "-"
    geometric_model(unit_cell, extr_dir, extr_dir_num_ele, extr_dir_cum_heights, extr_dir_recombine_ele)  
    #
    if !isempty(geom_export_paths)  # exporting UC model to the specified fomats
        for a_export_path in geom_export_paths
            verbose > 0 ? println("Saving geometric model tp $a_export_path") : "-"
            gmsh.write(a_export_path)
        end
    end
    #
    show_rve ? gmsh.fltk.run() : "-"
    #
    gmsh.finalize()
end


function make_lattice_cell_model(
    lattice_cell::AbstractLatticeCell;
    small_parameter::Float64=1e-06,
    export_paths::Union{String, Vector{String}}="None",  # joinpath(homedir(), "unit_cell.step"),
    mesh_periodicity::Bool=true,
    element_types::Tuple{Vararg{String}}=(),
    min_ele_size_factor::Float64=1.0,   # FIXME make default to take factor, instead of absolute values
    max_ele_size_factor::Float64=2.0,
    mesh_opt_algorithm::String="Netgen",
    node_renum_algorithm::String="RCM",  # give blank string "" for no renumbering
    show_mesh_stats::Bool=true,
    show_model::Bool=true,
    verbose::Int=1
)
    lc_dim::Int = 3
    xc, yc, zc = lattice_cell.centre
    lx, ly, lz = side_lengths(lattice_cell)
    small_side_length::Float64 = minimum([lx, ly, lz])
    lc_id = uppercase(lattice_cell.id)
    global ϵ = small_parameter
    #
    gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", 0)
    gmsh.option.set_number("General.Terminal", 0)
    gmsh.option.set_number("General.NoPopup", 0)
    gmsh.option.set_number("Geometry.OCCBoundsUseStl", 1)
    gmsh.model.add("unit_cell")
    # ============================================
    #               creating geometry
    # ============================================
    if isa(lattice_cell, SimpleCubicLattice)
        lattice_dim_tags = make_simple_cubic_lattice_geometry(
            xc, yc, zc, 
            lx, ly, lz, 
            lattice_cell.link_eqr_x, 
            lattice_cell.link_eqr_y, 
            lattice_cell.link_eqr_z
        )
    elseif lc_id == "FCC"
        @warn "FCC lattice not yet implemented"
    else
        throw(ArgumentError("Unknown lattice type: $(lattice_cell)"))
    end
    println(lattice_dim_tags)
    inclusions_pg_tags = gmsh.model.add_physical_group(lc_dim, [i[2] for i in lattice_dim_tags])
    gmsh.model.set_physical_name(lc_dim, inclusions_pg_tags, lc_id)
    # ============================================
    #               creating mesh
    # ============================================
    generate_mesh(
        lattice_cell,
        mesh_periodicity,
        element_types,
        min_ele_size_factor * small_side_length,
        max_ele_size_factor * small_side_length,
        mesh_opt_algorithm;
        verbose=verbose
    )
    show_mesh_stats && (verbose > 0) ? write_mesh_statistics() : "-"
    show_model ? gmsh.fltk.run() : "-"
    # ============================================
    #               exporting geometry model
    # ============================================
    if isa(export_paths, String)
        export_paths != "None" ? gmsh.write(export_paths) : nothing
    elseif isa(export_paths, Vector{String})
        for a_export_path in export_paths
            a_export_path != "None" ? gmsh.write(a_export_path) : nothing
        end
    end
    # ============================================
    #               exporting model mesh data
    # ============================================
    lc_geom_model_dim, lc_geom_model_tag = lattice_dim_tags[1]
    model_data::Dict{String,Any} = Dict{String, Any}()
    model_data["side_lengths"] = [lx, ly, lz]
    model_data["centre"] = [xc, yc, zc]
    model_data["mesh_data"] = get_mesh_data(
        lc_dim, 
        Dict("lattice" => [lc_geom_model_tag,],),
        node_renum_algorithm,
        verbose,
    )
    # 
    gmsh.finalize()
    # 
    return model_data
end

