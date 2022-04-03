
# elment Types
const CPE3 = "CPE3"
const CPE4 = "CPE4"
const CPS3 = "CPS3"
const CPS4 = "CPS4"
const C3D4 = "C3D4"
const C3D6 = "C3D6"
const C3D8 = "C3D8"


function apply_mesh_periodicity_constraints(
    uc::AbstractUnitCell
)
    dim::Int = dimension(uc)
    parents::Tuple = dim == 2 ? (:XLB, :YLB,) : (:XLB, :YLB, :ZLB)
    uc_side_lengths = side_lengths(uc)
    for a_parent in parents
        #
        affine_matrix = [
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
        ]
        #
        if a_parent == :XLB
            affine_matrix[4] = uc_side_lengths[1]
        elseif a_parent == :YLB
            affine_matrix[8] = uc_side_lengths[2]
        elseif a_parent == :ZLB
            affine_matrix[12] = uc_side_lengths[3]
        end
        #
        if isa(uc, UDC2D)
            pxl, pyl, pxu, pyu = buffer_bbox(uc, a_parent)
            pzl, pzu = 0.0, 0.0
        else
            pxl, pyl, pzl, pxu, pyu, pzu = buffer_bbox(uc, a_parent)
        end
        #
        parent_entities = gmsh.model.get_entities_in_bounding_box(pxl, pyl, pzl, pxu, pyu, pzu, dim-1)
        for (ap_dim, ap_tag) in parent_entities
            apxl, apyl, apzl, apxu, apyu, apzu = gmsh.model.get_bounding_box(ap_dim, ap_tag)
            apxl += (affine_matrix[4, ] - ϵ)
            apyl += (affine_matrix[8, ] - ϵ)
            apzl += (affine_matrix[12,] - ϵ)
            apxu += (affine_matrix[4, ] + ϵ)
            apyu += (affine_matrix[8, ] + ϵ)
            apzu += (affine_matrix[12,] + ϵ)
            chidren_entities = gmsh.model.get_entities_in_bounding_box(
                apxl, apyl, apzl, apxu, apyu, apzu, dim-1
            )
            for (ac_dim, ac_tag) in chidren_entities
                (acxl, acyl, aczl, acxu, acyu, aczu) = gmsh.model.get_bounding_box(ac_dim, ac_tag)
                if (
                    abs(acxl-apxl) < ϵ &&
                    abs(acyl-apyl) < ϵ &&
                    abs(aczl-apzl) < ϵ &&
                    abs(acxu-apxu) < ϵ &&
                    abs(acyu-apyu) < ϵ &&
                    abs(aczu-apzu) < ϵ
                )
                    gmsh.model.mesh.set_periodic(dim, [ac_tag,], [ap_tag,], affine_matrix)
                end
            end
        end
    end
end


function check_generated_ele_types(
    ele_types::Tuple,
    dim::Int64
)
    num_tet     = convert(Int64, gmsh.option.get_number("Mesh.NbTetrahedra"))  # C3D4
    num_hex     = convert(Int64, gmsh.option.get_number("Mesh.NbHexahedra"))  # C3D8
    num_prism   = convert(Int64, gmsh.option.get_number("Mesh.NbPrisms"))  # C3D6
    num_pyramids= convert(Int64, gmsh.option.get_number("Mesh.NbPyramids"))  #C3D5
    num_trihedra= convert(Int64, gmsh.option.get_number("Mesh.NbTrihedra")) 
    num_tri  = convert(Int64, gmsh.option.get_number("Mesh.NbTriangles"))  # CPS3 || CPE3
    num_quad = convert(Int64, gmsh.option.get_number("Mesh.NbQuadrangles"))  # CPS4 || CPE4
    if dim==3
        if !(:C3D4 in ele_types) && (num_tet > 0)
            @warn "C3D4 / Tet elements are generated which are not requested...!"
        end
        if !(:C3D6 in ele_types) && (num_prism > 0)
            @warn "C3D6 / Prism elements are generated which are not requested...!"
        end
        if !(:C3D8 in ele_types) && (num_hex > 0)
            @warn "C3D8 / Brick elements are generated which are not requested...!"
        end    
        if num_pyramids >0
            @warn "Pyramidal elements are generated which are not requested...!"
        end
        if num_trihedra >0
            @warn "Trihedra elements are generated which are not requested...!"
        end
    else
        if !(:CPS3 in ele_types || CPE3 in ele_types) && (num_tri > 0)
            @warn "Triangular elements are generated which are not requested...!"
        end
        if !(:CPS4 in ele_types || CPE4 in ele_types) && (num_quad > 0)
            @warn "Quad elements are generated which are not requested...!"
        end
    end
end


function get_mesh_statistics(dim::Int64=-1)::Dict{String, Int64}
    mesh_stats = Dict{String, Int64}()
    if (dim == -1 ) || (dim == 3)
        mesh_stats["num_tet"]     = convert(Int64, gmsh.option.get_number("Mesh.NbTetrahedra"))
        mesh_stats["num_hex"]     = convert(Int64, gmsh.option.get_number("Mesh.NbHexahedra"))
        mesh_stats["num_prism"]   = convert(Int64, gmsh.option.get_number("Mesh.NbPrisms"))
        mesh_stats["num_pyramids"]= convert(Int64, gmsh.option.get_number("Mesh.NbPyramids"))
        mesh_stats["num_trihedra"]= convert(Int64, gmsh.option.get_number("Mesh.NbTrihedra"))    
    end
    if (dim == -1 ) || (dim == 2)
        mesh_stats["num_tri"] = convert(Int64, gmsh.option.get_number("Mesh.NbTriangles"))
        mesh_stats["num_quad"] = convert(Int64, gmsh.option.get_number("Mesh.NbQuadrangles"))
    end
    #
    mesh_stats["total_nodes"] = convert(Int64, gmsh.option.get_number("Mesh.NbNodes"))
    #
    mesh_stats["num_3D_elements"] = (
        mesh_stats["num_tet"] +
        mesh_stats["num_hex"] +
        mesh_stats["num_prism"] +
        mesh_stats["num_pyramids"] +
        mesh_stats["num_trihedra"]
    )
    mesh_stats["num_2D_elements"] = (
        mesh_stats["num_tri"] +
        mesh_stats["num_quad"]
    )
    #
    return mesh_stats
end


function write_mesh_statistics(dim::Int64=-1)
    #
    mesh_stats = get_mesh_statistics(dim)
    #
    println(repeat("=", 50))
    printstyled("MESH STATISTICS\n"; color=:yellow, bold=true)
    println(repeat("=", 50))
    #
    println("--- NODES ---")
    println("Number of nodes ::                 \t", mesh_stats["total_nodes"])
    println(repeat("-", 10))
    println("--- 3D ELEMENTS ---")
    println("Number of Tetrahedra elements ::   \t", mesh_stats["num_tet"])
    println("Number of Hexahedral elements ::   \t", mesh_stats["num_hex"])
    println("Number of Prsim elements ::        \t", mesh_stats["num_prism"])
    println("Number of Pyramidal elements ::    \t", mesh_stats["num_pyramids"])
    println("Number of Trihedral elements ::    \t", mesh_stats["num_trihedra"])
    println(repeat("-", 20))
    println("TOTAL ::                           \t", mesh_stats["num_3D_elements"])
    println(repeat("-", 20))
    #

    println("--- 2D ELEMENTS ---")
    println("Number of Triangles ::             \t", mesh_stats["num_tri"])
    println("Number of Quadrangles ::           \t", mesh_stats["num_quad"])
    println(repeat("-", 20))
    println("TOTAL ::                           \t", mesh_stats["num_2D_elements"])
    println(repeat("-", 20))
    println(repeat("=", 50))
end

#  , , , mesh_min_size, mesh_max_size, mesh_opt_algorithm, show_mesh_stats

function generate_mesh(
    uc::AbstractUnitCell,
    mesh_periodicity::Bool,
    element_types::Tuple{Vararg{Symbol}},
    extr_dir_num_ele::Vector{Int64},
    min_ele_size::Float64,
    max_ele_size::Float64,
    mesh_opt_algorithm::String,
)
    #
    if mesh_periodicity
        apply_mesh_periodicity_constraints(uc)
    end
    #
    element_types = begin
        if (element_types[1:end] == (:DEFAULT, ) )
            if isa(uc, UDC2D)
                (CPS3,)
            elseif isa(uc, UDC3D)
                isempty(extr_dir_num_ele) ? (C3D4,) : (C3D6, C3D8,)
            elseif isa(uc, PRC)
                (C3D4,)
            end            
        else
            element_types
        end    
    end
    # 
    # Setting constraints to get the mesh of desired element types
    if any([i in (:CPS4, :CPE4, :C3D8) for i in element_types])
        gmsh.option.set_number("Mesh.RecombineAll", 1)
    end
    if (length(element_types)==1 && (element_types[1] == :CPE4 || element_types[1] == :CPS4 || element_types[1] == :C3D8))
        # Purely quad elements in the mesh
        gmsh.option.set_number("Mesh.Algorithm", 8)
        gmsh.option.set_number("Mesh.RecombinationAlgorithm", 3)
    end
    #
    gmsh.option.set_number("Mesh.MeshSizeMin", min_ele_size)
    gmsh.option.set_number("Mesh.MeshSizeMax", max_ele_size)
    #
    gmsh.model.mesh.generate(dimension(uc))
    gmsh.model.mesh.optimize(mesh_opt_algorithm, true)
    gmsh.model.mesh.remove_duplicate_nodes()
    gmsh.model.mesh.renumber_nodes()
    #
    check_generated_ele_types(element_types, dimension(uc))
    #
end


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
    geom_dim::Int,
    geometry_tags::Vector{Int},
)::Dict{String, Any}
    mesh_data = Dict{String, Any}()
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

    function _element_connectivity(gtags::Vector{Int64})::Dict{Int64, Matrix{Int64}}
        element_connectivity = Dict{Int64, Matrix{Int64}}()
        for ag_tag in gtags
            ele_types, ele_tags, ele_node_tags = gmsh.model.mesh.get_elements(geom_dim, ag_tag)
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

