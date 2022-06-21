

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


<<<<<<< HEAD
=======
# function write_element_database(
#     analysis_type::DataType,
#     all_node_set::Dict,
#     felset_data::Vector{FEPreProcessing.FiniteElementSet},
#     nodal_dataBase::Dict{Int,Matrix{Float64}},
#     strain_energies::Dict{Int,Float64},
# )::Dict{String,Matrix{Float64}}
#     #
#     Element_DataBase::Dict{String,Matrix{Float64}} = Dict{String,Matrix{Float64}}()
#     #
#     num_variables = length(variables)
#     #
#     numFele::Dict = FEPreProcessing.get_num_ele(felset_data)
#     ele_counters::Dict = Dict{String,Int}()
#     for (a_fele_type, num_felements) in numFele
#         num_nodes_per_ele = FEPreProcessing.get_num_nodes_per_ele(a_fele_type)
#         a_ele_type = split(string(a_fele_type), ".")[end]
#         Element_DataBase[a_ele_type] = zeros(Float64, (1 + (num_variables * num_nodes_per_ele)), num_felements)
#         ele_counters["max_"*a_ele_type] = num_felements
#         ele_counters[a_ele_type] = 0
#     end

#     for a_fele_set in felset_data
#         # running over each material/finte element set
#         felements = a_fele_set.elements
#         for a_fele in felements
#             # running over each finite element of the current finite element set
#             a_ele_tag = a_fele.tag
#             node_tags = a_fele.node_tags
#             num_nodes_per_ele = length(node_tags)
#             a_ele_data = zeros(Float64, num_variables * num_nodes_per_ele + 1)
#             for (node_idx, a_nodeTag) in enumerate(node_tags)
#                 # Running over each node of the current element
#                 a_node_x, a_node_y, a_node_z = all_node_set[a_nodeTag]
#                 a_node_field_data = nodal_dataBase[a_nodeTag]  # 15x1 matrix
#                 #
#                 # ---x---y---z---ux---uy---uz---E11---E22---E33---E23---E31---E12---S11---S22---S33---S23---S31---S12---SE
#                 # ==============================
#                 #  adding nodal coordinate data
#                 # ==============================
#                 a_ele_data[node_idx] = a_node_x
#                 a_ele_data[node_idx+num_nodes_per_ele] = a_node_y
#                 a_ele_data[node_idx+(2*num_nodes_per_ele)] = a_node_z
#                 # ==============================
#                 #  adding nodal deformation data
#                 # ==============================
#                 a_ele_data[node_idx+(3*num_nodes_per_ele)] = a_node_field_data[1, 1]
#                 a_ele_data[node_idx+(4*num_nodes_per_ele)] = a_node_field_data[2, 1]
#                 a_ele_data[node_idx+(5*num_nodes_per_ele)] = a_node_field_data[3, 1]
#                 # ==============================
#                 #  adding nodal strain data
#                 # ==============================
#                 a_ele_data[node_idx+(6*num_nodes_per_ele)] = a_node_field_data[4, 1]
#                 a_ele_data[node_idx+(7*num_nodes_per_ele)] = a_node_field_data[5, 1]
#                 a_ele_data[node_idx+(8*num_nodes_per_ele)] = a_node_field_data[6, 1]
#                 a_ele_data[node_idx+(9*num_nodes_per_ele)] = a_node_field_data[7, 1]
#                 a_ele_data[node_idx+(10*num_nodes_per_ele)] = a_node_field_data[8, 1]
#                 a_ele_data[node_idx+(11*num_nodes_per_ele)] = a_node_field_data[9, 1]
#                 # ==============================
#                 #  adding nodal Stress data
#                 # ==============================
#                 a_ele_data[node_idx+(12*num_nodes_per_ele)] = a_node_field_data[10, 1]
#                 a_ele_data[node_idx+(13*num_nodes_per_ele)] = a_node_field_data[11, 1]
#                 a_ele_data[node_idx+(14*num_nodes_per_ele)] = a_node_field_data[12, 1]
#                 a_ele_data[node_idx+(15*num_nodes_per_ele)] = a_node_field_data[13, 1]
#                 a_ele_data[node_idx+(16*num_nodes_per_ele)] = a_node_field_data[14, 1]
#                 a_ele_data[node_idx+(17*num_nodes_per_ele)] = a_node_field_data[15, 1]
#                 # Von-Mises Stress
#                 a_ele_data[node_idx+(18*num_nodes_per_ele)] = get_VonMises_stress(a_node_field_data[10:15, 1])
#             end
#             # ==============================
#             #  adding element strain energy
#             # ==============================
#             a_ele_data[end] = strain_energies[a_ele_tag]
#             a_ele_data = reshape(a_ele_data, :, 1)
#             #
#             # ==============================
#             #  update Element_DataBase
#             # ==============================
#             a_ele_type = split(string(typeof(a_fele)), ".")[end]
#             ele_counters[a_ele_type] += 1
#             ele_count = ele_counters[a_ele_type]
#             Element_DataBase[a_ele_type][:, ele_count] = a_ele_data
#         end
#     end
#     return Element_DataBase
# end

# function plot_a_field(
#     field_tag::String,
#     field_data::Dict{String,Matrix{Float64}},
#     plot_path::String;
#     field_type::String="S"
# )
#     gmsh.initialize()
#     field_view = gmsh.view.add(field_tag)
#     for (a_ele_type, ele_data) in field_data
#         FieldFele_Type = if a_ele_type in ("C3D8", "H", 5)
#             field_type * "H"
#         elseif a_ele_type in ("C3D86", "I", 6)
#             field_type * "I"
#         end
#         #
#         gmsh.view.add_list_data(
#             field_view,
#             FieldFele_Type,
#             size(ele_data, 2),
#             reshape(ele_data, :),
#         )
#     end
#     # gmsh.fltk.initialize()
#     # gmsh.write(plot_path)
#     gmsh.view.write(field_view, plot_path)
#     gmsh.view.remove(field_view)
#     gmsh.clear()
#     gmsh.finalize()
# end

# function plot_elastic_fields(
#     field_data::Dict{String,Matrix{Float64}},
#     plot_dir::String;
#     field_type::String="S",
#     plot_fname_suffix::String=""
# )
#     fields = [
#         "UX", "UY", "UZ",
#         "EPSILON-11", "EPSILON-22",
#         "EPSILON-33", "EPSILON-23",
#         "EPSILON-31", "EPSILON-12",
#         "SIGMA-11", "SIGMA-22",
#         "SIGMA-33", "SIGMA-23",
#         "SIGMA-31", "SIGMA-12",
#         "Von-Mises-Stress",
#     ]
#     for (idx, a_field) in enumerate(fields)
#         a_field_data::Dict{String,Matrix{Float64}} = Dict{String,Matrix{Float64}}()
#         for (aele_type, aele_type_data) in field_data
#             num_npe = if aele_type == "C3D8"  # num_nodes_per_ele
#                 8
#             elseif aele_type == "C3D6"
#                 6
#             end
#             #
#             coor_slicer = 1:(3*num_npe)
#             field_data_slicer = (1+(num_npe*(2+idx))):(num_npe*(3+idx))
#             a_field_data[aele_type] = aele_type_data[[coor_slicer; field_data_slicer], :]
#         end
#         plot_a_field(
#             a_field,
#             a_field_data,
#             joinpath(plot_dir, plot_fname_suffix * "_" * a_field * "_.pos");
#             field_type=field_type
#         )
#     end
# end


# function prepare_gmshpp_ele_data(
#     uc_odb::Dict{String,Any},
#     load_case_ID::String,
#     nodal_fvariables::Dict{String,Int},
#     element_fvariables::Dict{String,Int},
# )
#     nodal_coor = uc_odb["nodal_coordinates"]
#     felement_conn = uc_odb["element_connectivity"]
#     field_data = uc_odb[load_case_ID]
#     node_field_data = field_data["node_odb"]
#     element_field_data = field_data["element_odb"]
#     num_nodal_fvar = length(nodal_fvar_indices)
#     num_ele_fvar = length(element_fvar_indices)
#     #
#     felement_data_matrix::Dict{Int,Matrix{Float64}} = Dict{Int,Matrix{Float64}}()
#     for (a_felt, a_felt_conn) in felement_conn  # felt => finite element type
#         num_nodes_per_fele = size(a_felt_conn, 1) - 1
#         num_felements = size(a_felt_conn, 2)
#         a_felt_data_matrix = zeros(Float64, 3 + num_nodal_fvar + num_ele_fvar, num_felements)
#         # run over each element of the current type
#         for i in 1:num_felements
#             aelem_tag = a_felt_conn[1, i]
#             ntags = a_felt_conn[2:end, i]
#             #
#             # collecting coordinate and field data information for the current element
#             aele_node_coor = zeros(Float64, num_nodes_per_fele, 3)
#             aele_node_field_data = zeros(Float64, num_nodes_per_fele, num_nodal_fvar)
#             for (anode_idx, anode_tag) in enumerate(ntags)
#                 aele_node_coor[anode_idx, :] = nodal_coor[anode_tag]
#                 for an_fvar_idx in values(nodal_fvariables)
#                     aele_node_field_data[anode_idx, an_fvar_idx] = node_field_data[anode_tag][an_fvar_idx]
#                 end
#             end
#             aele_node_coor = reshape(aele_node_coor, :)
#             aele_node_field_data = reshape(aele_node_field_data, :)
#             a_felt_data_matrix[1:(3+num_nodal_fvar), i] = [aele_node_coor; aele_node_field_data]
#             #
#             for (k, aele_fvar_idx) in enumerate(values(element_fvariables))
#                 a_felt_data_matrix[(3+num_nodal_fvar+k), i] = element_field_data[aelem_tag][aele_fvar_idx]
#             end
#             #
#         end
#         felement_data_matrix[a_felt] = a_felt_data_matrix
#     end
# end
>>>>>>> 3cf2f09c55b4c940a9b6dbfd7467e9ef3149dce3
