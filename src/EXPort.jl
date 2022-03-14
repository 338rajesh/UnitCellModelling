

"""
    get_mesh_data()

Returns a dictionary of mesh data with the following key-values pairs

- `"all_node_tags"` => `Vector{Int}`
- `"all_node_coordinates"` => `Matrix{Float64}`
- `"matrix_element_connectivity"` => `Matrix{Int}`
- `"inclusions_element_connectivity"` => `Matrix{Int}`
- `"total_nodes"` => Int
- `"num_2D_elements"` => Int
- `"num_3D_elements"` => Int

"""
function get_mesh_data(
    uc::AbstractUnitCell,
    geometry_tags::Vector{Int},
)::Dict{String, Any}
    mesh_data = Dict{String, Any}()
    uc_dim::Int = dimension(uc)
    # -----------------------------------
    #   Nodal data collection
    # -----------------------------------
    nt, nc, _ = gmsh.model.mesh.get_nodes()
    mesh_data["all_node_tags"] = convert(Vector{Int64}, nt)  # Vector{Int64}
    mesh_data["all_node_coordinates"] = reshape(nc, 3, :)  # Matrix{Float64}
    #
    # -----------------------------------
    #   Collecting element connectivity
    # -----------------------------------
    element_connectivity::Dict{Int64, Matrix{Int64}} = Dict()
    ele_tags::Vector{Vector{Int}} = Vector{Int}[]
    ele_node_tags::Vector{Vector{Int}} = Vector{Int}[]
    ele_types::Vector{Int} = Int[]
    elt_num_nodes::Int = 0

    function _element_connectivity(gtags::Vector{Int64})
        element_connectivity = Dict{Int64, Matrix{Int64}}()
        for ag_tag in gtags
            ele_types, ele_tags, ele_node_tags = gmsh.model.mesh.get_elements(uc_dim, ag_tag)
            #
            for (etk, ael_type) in enumerate(ele_types)
                _, _, _, elt_num_nodes, _, _ = gmsh.model.mesh.get_element_properties(ael_type)
                elt_ele_connectivity = [
                    reshape(ele_tags[etk], 1, :);
                    reshape(ele_node_tags[etk], elt_num_nodes, :)
                ]
                # if element type (ael_type,) connectivity data is exisiting then append new data.
                if ael_type in keys(element_connectivity)
                    elt_ele_connectivity = hcat(element_connectivity[ael_type], elt_ele_connectivity)
                end
                #
                element_connectivity[ael_type] = elt_ele_connectivity                
            end    
        end
        return element_connectivity
    end
    mesh_data["matrix_element_connectivity"] = _element_connectivity(geometry_tags[1:1])
    mesh_data["inclusions_element_connectivity"] = _element_connectivity(geometry_tags[2:end])
    #
    # Adding mesh Statistics like number of nodes, elements....etc.
    merge!(mesh_data, get_mesh_statistics())
    #
    return mesh_data
end


function export_geometry(
    options::Dict
)
    if :geom_export_paths in keys(options)
        for a_export_path in options[:geom_export_paths]
            gmsh.write(a_export_path)
        end
    end
end

