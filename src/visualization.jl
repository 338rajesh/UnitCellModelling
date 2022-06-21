

function plot_nodal_fields(
    uc_odb::Dict{String,Any},
    load_case_IDs::Vector{String},
    nodal_fvariables::Dict{String,Int};
    plots_dir::String=mkpath(joinpath(homedir(), "uc_plots")),
    field_type::String="S",
    plot_type::String=".png",
)
    nodal_coor = uc_odb["nodal_coordinates"]
    felement_conn = uc_odb["element_connectivity"]
    #
    gmsh.initialize()
    for a_load_case_ID in load_case_IDs
        node_field_data = uc_odb[a_load_case_ID]["node_odb"]
        for (nodal_fvar_name, nodal_fvar_index) in nodal_fvariables
            println("Plotting $nodal_fvar_name field for $a_load_case_ID load case.")
            field_view = gmsh.view.add(a_load_case_ID * ": " * nodal_fvar_name)
            for (a_felt, a_felt_conn) in felement_conn  # felt => finite element type
                npe = size(a_felt_conn, 1) - 1  # number of nodes per element
                num_felements = size(a_felt_conn, 2)
                felement_type_code = if a_felt in ("C3D8", "H", 5)
                    "H"
                elseif a_felt in ("C3D86", "I", 6)
                    "I"
                end
                # run over each element of the current type
                fvar_data = zeros(Float64, 4 * num_felements * npe)
                for i in 1:num_felements
                    # collecting coordinate and field data information for the current element
                    afele_start_idx = 4 * (i-1) * npe
                    for j in 1:npe
                        ntag = a_felt_conn[1+j, i]
                        # extract node coordinates
                        ax, ay, az = nodal_coor[ntag]
                        # extract node field variable
                        anode_fvar = node_field_data[ntag][nodal_fvar_index]
                        fvar_data[afele_start_idx + 0*npe + j] = ax
                        fvar_data[afele_start_idx + 1*npe + j] = ay
                        fvar_data[afele_start_idx + 2*npe + j] = az
                        fvar_data[afele_start_idx + 3*npe + j] = anode_fvar
                    end
                end
                #
                gmsh.view.add_list_data(
                    field_view,
                    field_type * felement_type_code,
                    num_felements,
                    fvar_data,
                )
            end
            plot_path = joinpath(plots_dir, a_load_case_ID * "_" * nodal_fvar_name * plot_type)
            if plot_type == ".pos"
                gmsh.view.write(field_view, plot_path)
            else
                gmsh.fltk.initialize()
                gmsh.write(plot_path)
            end
            gmsh.view.remove(field_view)
            gmsh.clear()
        end
    end
    gmsh.finalize()        
end


