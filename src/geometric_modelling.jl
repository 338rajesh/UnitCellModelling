

function add_bbox(
    unit_cell::Union{UDC2D,UDC3D}
)::Int
    if isa(unit_cell, UDC2D)
        xlb, ylb, xub, yub = unit_cell.bbox
        bbox_tag = gmsh.model.occ.add_rectangle(xlb, ylb, 0.0, xub - xlb, yub - ylb,)
    else
        xlb, ylb, zlb, xub, yub, _ = unit_cell.bbox
        bbox_tag = gmsh.model.occ.add_rectangle(xlb, ylb, zlb, xub - xlb, yub - ylb,)
    end
    return bbox_tag
end


function add_inclusions(
    unit_cell::Union{UDC2D,UDC3D}
)::Vector{Int}
    uc_inclusions::Dict{String,Matrix{Float64}} = unit_cell.inclusions
    z_min = isa(unit_cell, UDC2D) ? 0.0 : unit_cell.bbox.zlb
    incl_tags::Vector{Int} = Int[]
    for (incl_shape_id, incl_data) in uc_inclusions
        incl_tags = [
            incl_tags
            begin
                if uppercase(incl_shape_id) == "CIRCLE"
                    add_circular_discs(incl_data, z_min)
                elseif uppercase(incl_shape_id) == "ELLIPSE"
                    add_elliptical_discs(incl_data, z_min)
                elseif uppercase(incl_shape_id) == "CAPSULE"
                    add_capsular_discs(incl_data, z_min)
                elseif uppercase(incl_shape_id) == "RECTANGLE"
                    add_rectangular_discs(incl_data, z_min)
                elseif uppercase(incl_shape_id) == "RPOLYGON"
                    add_rpolygons(incl_data, z_min)
                end
            end
        ]
    end
    return incl_tags
end


"""
    get_matrix_phase_tag()

    with positive buffer(pb) and negative buffer(nb), surfaces
    belonging to a slightly enlarged and slightly shrinked bboxes than
    that of the unit cell.

"""
function get_matrix_phase_tag(
    uc::AbstractUnitCell,
    entity_dim::Int,
)::Int
    if isa(uc, UDC2D)
        xl, yl, xu, yu = uc.bbox
        zl, zu = 0.0, 0.0
    elseif isa(uc, UDC3D)
        if entity_dim == 2
            xl, yl, zl, xu, yu, _ = uc.bbox
            zu = zl
        else
            xl, yl, zl, xu, yu, zu = uc.bbox        
        end
    elseif isa(uc, PRC)
        xl, yl, zl, xu, yu, zu = uc.bbox
    end
    #
    res = Int[]
    entities = gmsh.model.get_entities_in_bounding_box(
        xl-ϵ, yl-ϵ, zl-ϵ, xu+ϵ, yu+ϵ, zu+ϵ, entity_dim
    )
    for (ae_dim, ae_tag) in entities
        ae_xl, ae_yl, ae_zl, ae_xu, ae_yu, ae_zu = gmsh.model.get_bounding_box(ae_dim, ae_tag)
        if (
            abs(xl - ae_xl) < ϵ && abs(yl - ae_yl) < ϵ && abs(zl - ae_zl) < ϵ &&
            abs(xu - ae_xu) < ϵ && abs(yu - ae_yu) < ϵ && abs(zu - ae_zu) < ϵ
        )
            push!(res, ae_tag)
        end
    end
    if length(res) == 1
        return res[1]
    elseif length(res) == 0
        @error "Unable to find an entity with the given bbox"
        exit(1)
    elseif length(res) > 1
        @error "More than one entity is found with the given bbox"
        exit(1)
    end
end


function make_udc_transverse_section(
    uc::UDC,
)::Vector{Int}
    # ===========================================
    #           creating raw 2D parts
    # ===========================================
    raw_matrix_tag = add_bbox(uc)
    raw_inclusions_tags = add_inclusions(uc)
    gmsh.model.occ.synchronize()
    # ===========================================
    #           BOOLEAN OPERATIONS
    # ===========================================
    fragments, _ = gmsh.model.occ.fragment(
        (2, raw_matrix_tag),
        [(2, i) for i in raw_inclusions_tags]
    )
    gmsh.model.occ.synchronize() 
    
    if isa(uc, UDC2D)
        xlb_buff, ylb_buff, xub_buff, yub_buff = buffer_bbox(uc, :ALL, ϵ)
        zlb_buff, zub_buff = 0.0, 0.0
    else
        xlb_buff, ylb_buff, zlb_buff, xub_buff, yub_buff, zub_buff = buffer_bbox(uc, :ALL, ϵ)
    end

    # Checking and removing, if there are any surfaces outside unit cell.
    surfaces_in_uc_bounds = gmsh.model.get_entities_in_bounding_box(
        xlb_buff, ylb_buff, zlb_buff, xub_buff, yub_buff, zub_buff, 2
    )
    outer_surface_fragments = [i for i in fragments if !(i in surfaces_in_uc_bounds)]
    gmsh.model.occ.remove(outer_surface_fragments, true)
    gmsh.model.occ.synchronize()

    # Checking and removing, if there are any curves outside unit cell.
    all_curves = gmsh.model.get_entities(1)
    curves_within_uc_bounds = gmsh.model.get_entities_in_bounding_box(
        xlb_buff, ylb_buff, zlb_buff, xub_buff, yub_buff, zub_buff, 1
    )
    exterior_curves = [a_curve for a_curve in all_curves if !(a_curve in curves_within_uc_bounds)]
    gmsh.model.occ.remove(exterior_curves, true)
    gmsh.model.occ.synchronize()

    # Checking and removing, if there are any points outside unit cell.
    all_points = gmsh.model.get_entities(0)
    # points_within_uc_bounds = gmsh.model.get_entities_in_bounding_box(
    #     xlb_buff, ylb_buff, zlb_buff, xub_buff, yub_buff, zub_buff, 0
    # )
    # exterior_points = [a_pnt for a_pnt in all_points if !(a_pnt in points_within_uc_bounds)]

    # println("Number of exterior_points ", length(exterior_points))
    gmsh.model.occ.remove(all_points)
    gmsh.model.occ.synchronize()

    # ===========================================
    #           FINDING THE MATRIX TAG
    # ===========================================
    #

    matrix_tag = get_matrix_phase_tag(uc, 2)

    inclusion_tags = [i for (_, i) in surfaces_in_uc_bounds if i != matrix_tag]

    return [matrix_tag; inclusion_tags]
end


function geometric_model(
    uc::UDC,
    extr_dir::String,
    extr_dir_num_ele::Vector{Int64},
    extr_dir_cum_heights::Vector{Float64},
    extr_dir_recombine_ele::Bool,
)::Vector{Int}
    uc_cross_section_surfaces = make_udc_transverse_section(uc)
    if isa(uc, UDC2D)
        return uc_cross_section_surfaces
    else
        if extr_dir == "XY->Z"
            dx, dy, dz = 0.0, 0.0, uc.bbox.zub - uc.bbox.zlb
        elseif extr_dir == "ZX->Y"
            dx, dy, dz = 0.0, uc.bbox.yub - uc.bbox.ylb, 0.0
        elseif extr_dir == "YZ->X"
            dx, dy, dz = uc.bbox.xub - uc.bbox.xlb, 0.0, 0.0
        end
        extrude_dimTags = begin
            if !isempty(extr_dir_num_ele)
                gmsh.model.occ.extrude(
                    [(2, i) for i in uc_cross_section_surfaces],
                    dx, dy, dz,
                    extr_dir_num_ele,
                    extr_dir_cum_heights,
                    extr_dir_recombine_ele,
                )
            else
                gmsh.model.occ.extrude(
                    [(2, i) for i in uc_cross_section_surfaces],
                    dx, dy, dz,
                )
            end
        end
        gmsh.model.occ.synchronize()
        unit_cell_vol_tags = [i_tag for (i_dim, i_tag) in extrude_dimTags if i_dim == 3]
        matrix_tag = get_matrix_phase_tag(uc, 3)
        inclusion_tags = [i for i in unit_cell_vol_tags if i != matrix_tag]
        return [matrix_tag; inclusion_tags]
    end
end

