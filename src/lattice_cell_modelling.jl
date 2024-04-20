function trim_to_bbox(
    obj_dim_tags::Vector{Tuple{Int32, Int32}},
    bbox_dim_tags::Vector{Tuple{Int32, Int32}},
)::Vector{Tuple{Int64, Int64}}
    # Preparing the tool
    b_box_2_dimTags = gmsh.model.occ.copy(bbox_dim_tags)
    xc, yc, zc = gmsh.model.occ.get_center_of_mass(b_box_2_dimTags[1][1], b_box_2_dimTags[1][2])
    gmsh.model.occ.dilate(b_box_2_dimTags, xc, yc, zc, 2.0, 2.0, 2.0)
    trimming_tool_dim_tags, _ = gmsh.model.occ.cut(b_box_2_dimTags, bbox_dim_tags)
    gmsh.model.occ.synchronize()
    # Cut the object with the tool
    obj_dim_tags, _ = gmsh.model.occ.cut(obj_dim_tags, trimming_tool_dim_tags)
    gmsh.model.occ.synchronize()
    return obj_dim_tags
end



"""
    ler::Float64, Link Equivalent Radius
    xc, yc, zc::Float64, Center of the lattice
    lx, ly, lz::Float64, Lattice bounding box dimensions
"""
function make_simple_cubic_lattice_geometry(
    xc::Float64, yc::Float64, zc::Float64,
    lx::Float64, ly::Float64, lz::Float64,
    eqr_x::Float64, eqr_y::Float64, eqr_z::Float64
)
    # Adding bounding box
    b_box = gmsh.model.occ.add_box(0.0, 0.0, 0.0, lx, ly, lz)
    # Adding cylinders Parallel to X-Axis
    p12 = gmsh.model.occ.add_cylinder(0.0, 0.0, 0.0, lx, 0.0, 0.0, eqr_x)
    p43 = gmsh.model.occ.add_cylinder(0.0, 0.0, lz, lx, 0.0, 0.0, eqr_x)
    p87 = gmsh.model.occ.add_cylinder(0.0, ly, lz, lx, 0.0, 0.0, eqr_x)
    p56 = gmsh.model.occ.add_cylinder(0.0, ly, 0.0, lx, 0.0, 0.0, eqr_x)
    # Adding cylinders Parallel to Y-Axis
    p15 = gmsh.model.occ.add_cylinder(0.0, 0.0, 0.0, 0.0, ly, 0.0, eqr_y)
    p26 = gmsh.model.occ.add_cylinder(lx, 0.0, 0.0, 0.0, ly, 0.0, eqr_y)
    p37 = gmsh.model.occ.add_cylinder(lx, 0.0, lz, 0.0, ly, 0.0, eqr_y)
    p48 = gmsh.model.occ.add_cylinder(0.0, 0.0, lz, 0.0, ly, 0.0, eqr_y)
    # Adding cylinders Parallel to Z-Axis
    p14 = gmsh.model.occ.add_cylinder(0.0, 0.0, 0.0, 0.0, 0.0, lz, eqr_z)
    p23 = gmsh.model.occ.add_cylinder(lx, 0.0, 0.0, 0.0, 0.0, lz, eqr_z)
    p67 = gmsh.model.occ.add_cylinder(lx, ly, 0.0, 0.0, 0.0, lz, eqr_z)
    p58 = gmsh.model.occ.add_cylinder(0.0, ly, 0.0, 0.0, 0.0, lz, eqr_z)
    gmsh.model.occ.synchronize()
    # Fusing cylinders
    lattice_dimTags, _ = gmsh.model.occ.fuse(
        [[3, i] for i in [p12, p43, p87, p56, ]],
        [[3, i] for i in [p15, p26, p37, p48, p14, p23, p67, p58, ]],
    )
    gmsh.model.occ.synchronize()
    # 
    lattice_dimTags = trim_to_bbox(lattice_dimTags, [(Int32(3), b_box)])
    #
    gmsh.model.occ.translate(lattice_dimTags, xc-0.5*lx, yc-0.5*ly, zc-0.5*lz)
    gmsh.model.occ.synchronize()
    return lattice_dimTags
end


