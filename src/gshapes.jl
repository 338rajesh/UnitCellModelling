# ==========================================
#          REGULAR 2D SHAPES
# ==========================================


function make_rpolygon(
    xyz_centre::NTuple{3,Float64},
    theta::Float64,
    side_len::Float64,
    rf::Float64,
    num_sides::Real;
    normal_vec::NTuple{3,Float64}=(0.0, 0.0, 1.0)
)::Int
    num_sides = convert(Int64, num_sides)
    xc, yc, zc = xyz_centre
    ax, ay, az = normal_vec
    α::Float64 = π / num_sides
    circum_radius::Float64 = side_len / (2.0 * sin(α))
    h::Float64 = rf / cos(α)
    # 
    pointsTags = []
    for asn in 1:num_sides
        θ_i = 2.0 * α * (asn - 1.0)
        o_dimTag = gmsh.model.occ.add_point(xc + circum_radius - h, yc, zc)
        p_dimTag = gmsh.model.occ.add_point(xc + circum_radius - h + rf * cos(α), yc + (rf * sin(α)), zc)
        q_dimTag = gmsh.model.occ.add_point(xc + circum_radius - h + rf * cos(α), yc - (rf * sin(α)), zc)
        gmsh.model.occ.rotate([(0, o_dimTag)], xc, yc, zc, ax, ay, az, θ_i)
        gmsh.model.occ.rotate([(0, p_dimTag)], xc, yc, zc, ax, ay, az, θ_i)
        gmsh.model.occ.rotate([(0, q_dimTag)], xc, yc, zc, ax, ay, az, θ_i)
        push!(pointsTags, (p_dimTag, o_dimTag, q_dimTag))
    end
    curvesTags = []
    for asn in 1:num_sides
        arcTag = gmsh.model.occ.add_circle_arc(pointsTags[asn]...)
        if asn <= num_sides - 1
            lineTag = gmsh.model.occ.add_line(pointsTags[asn][1], pointsTags[asn+1][3])
        elseif asn == num_sides
            lineTag = gmsh.model.occ.add_line(pointsTags[asn][1], pointsTags[1][3])
        end
        curvesTags = [curvesTags; [arcTag, lineTag]]
    end
    rpolygon_wireTag = gmsh.model.occ.add_curve_loop(curvesTags)
    rpolygon_discTag = gmsh.model.occ.add_plane_surface([rpolygon_wireTag])
    gmsh.model.occ.rotate((2, rpolygon_discTag), xc, yc, zc, ax, ay, az, theta)
    # gmsh.model.occ.synchronize()
    #
    return rpolygon_discTag
end


function make_capsule(
    xyz::NTuple{3,Float64},
    theta::Float64,
    smj_len::Float64,
    smn_len::Float64,
)::Int
    xc, yc, zc = xyz
    p1 = gmsh.model.occ.add_point(xc + smj_len - smn_len, yc + smn_len, zc)
    p2 = gmsh.model.occ.add_point(xc - smj_len + smn_len, yc + smn_len, zc)
    p3 = gmsh.model.occ.add_point(xc - smj_len + smn_len, yc - smn_len, zc)
    p4 = gmsh.model.occ.add_point(xc + smj_len - smn_len, yc - smn_len, zc)
    c1 = gmsh.model.occ.add_point(xc + smj_len - smn_len, yc, zc)
    c2 = gmsh.model.occ.add_point(xc - smj_len + smn_len, yc, zc)
    # 
    line_1 = gmsh.model.occ.add_line(p1, p2)
    arc_1 = gmsh.model.occ.add_circle_arc(p3, c2, p2)
    line_2 = gmsh.model.occ.add_line(p3, p4)
    arc_2 = gmsh.model.occ.add_circle_arc(p1, c1, p4)
    # 
    caps_wire_tag = gmsh.model.occ.add_curve_loop([line_1, arc_1, line_2, arc_2])
    caps_disc_tag = gmsh.model.occ.add_plane_surface([caps_wire_tag])
    gmsh.model.occ.rotate((2, caps_disc_tag), xc, yc, zc, 0.0, 0.0, 1.0, theta)
    return caps_disc_tag
end


function make_ellipse(
    xyz::NTuple{3,Float64},
    theta::Float64,
    smj_len::Float64,
    smn_len::Float64,
)::Int
    xc, yc, zc = xyz
    ell_curve_tag = gmsh.model.occ.add_ellipse(xc, yc, zc, smj_len, smn_len)
    ell_wire_tag = gmsh.model.occ.add_wire([ell_curve_tag])
    ell_disc_tag = gmsh.model.occ.add_plane_surface([ell_wire_tag])
    gmsh.model.occ.rotate((2, ell_disc_tag), xc, yc, zc, 0.0, 0.0, 1.0, theta)
    return ell_disc_tag
end


function make_rectangle(
    xyz::NTuple{3,Float64},
    theta::Float64,
    smj_len::Float64,
    smn_len::Float64,
    corner_radius::Float64
)::Int
    xc, yc, zc = xyz
    rect_tag = gmsh.model.occ.add_rectangle(xc - smj_len, yc - smn_len, zc, 2.0 * smj_len, 2.0 * smn_len, -1, corner_radius)
    gmsh.model.occ.rotate((2, rect_tag), xc, yc, zc, 0.0, 0.0, 1.0, theta)
    return rect_tag
end


function make_cshape(
    xyz_centre::NTuple{3,Float64},
    θ1::Float64,  # Orientation of starting tip.
    ro::Float64,  # outer radius
    ri::Float64,  # inner radius
    θc::Float64,  # included angle  # FIXME check for alpha greater than π
)::Int
    xc, yc, zc = xyz_centre
    r_mean::Float64 = 0.5 * (ro + ri)
    r_tip::Float64 = 0.5 * (ro - ri)
    θ2::Float64 = θ1 + (0.5*θc)
    θ3::Float64 = θ1 + (1.0*θc)
    #
    add_a_point(
        x_::Float64, y_::Float64, z_::Float64, tht::Float64, r_::Float64
    )::Int = gmsh.model.occ.add_point(x_ + (r_ * cos(tht)), y_ + (r_ * sin(tht)), z_)
    
    c0::Int = add_a_point(xc, yc, zc, 0.0, 0.0)
    c4::Int = add_a_point(xc, yc, zc, θ3, r_mean)
    c8::Int = add_a_point(xc, yc, zc, θ1, r_mean)
    #
    # Adding points
    p_1::Int = add_a_point(xc, yc, zc, θ1, ro)
    p_2::Int = add_a_point(xc, yc, zc, θ2, ro)
    p_3::Int = add_a_point(xc, yc, zc, θ3, ro)
    p_4::Int = add_a_point(xc+(r_mean*cos(θ3)), yc+(r_mean*sin(θ3)), zc, θ3+(0.5*π), r_tip)  # tip-1
    p_5::Int = add_a_point(xc, yc, zc, θ3, ri)
    p_6::Int = add_a_point(xc, yc, zc, θ2, ri)
    p_7::Int = add_a_point(xc, yc, zc, θ1, ri)
    p_8::Int = add_a_point(xc+(r_mean*cos(θ1)), yc+(r_mean*sin(θ1)), zc, θ1+(1.5*π), r_tip)  # tip-2
    #
    gmsh.model.occ.synchronize()
    # Adding circular arcs and making wire
    cshape_wire_tag = gmsh.model.occ.add_curve_loop([
        gmsh.model.occ.add_circle_arc(p_1, c0, p_2),
        gmsh.model.occ.add_circle_arc(p_2, c0, p_3),
        gmsh.model.occ.add_circle_arc(p_3, c4, p_4),
        gmsh.model.occ.add_circle_arc(p_4, c4, p_5),
        gmsh.model.occ.add_circle_arc(p_5, c0, p_6),
        gmsh.model.occ.add_circle_arc(p_6, c0, p_7),
        gmsh.model.occ.add_circle_arc(p_7, c8, p_8),
        gmsh.model.occ.add_circle_arc(p_8, c8, p_1),
    ])
    #
    # Making Surface
    cshape_disc_tag = gmsh.model.occ.add_plane_surface([cshape_wire_tag,])
    return cshape_disc_tag
end


function make_nlobeshape(
    xyz_centre::NTuple{3,Float64},  # centre of the n-lobe
    theta::Float64,  # Orientation of the initial/reference tip
    ro::Float64,  # outer radius
    rl::Float64,  # lobe radius
    num_lobes::Int;
    normal_vec::NTuple{3, Float64}=(0.0, 0.0, 1.0),
)::Int
    α::Float64 = π / num_lobes
    θ::Float64 = asin(sin(α) * ((ro-rl)/(2.0*rl)))
    # ======================
    #   MAKING UNIT TIP
    # ======================
    # Evaluating points
    C = xyz_centre
    C1 = C .+ (ro-rl, 0.0, 0.0)
    C2 = C1 .+ (2.0*rl*cos(α+θ), 2.0*rl*sin(α+θ), 0.0)
    P1 = C1 .+ (rl, 0.0, 0.0)
    P2 = C1 .+ (rl*cos(α+θ), rl*sin(α+θ), 0.0)
    P3 = C2 .+ (-rl*cos(α), -rl*sin(α), 0.0)
    # Making Gmsh points
    add_nl_points(tpl) = gmsh.model.occ.add_point(tpl[1], tpl[2], tpl[3])
    # gm_C::Int    = add_nl_points(C)
    gm_C1::Int    = add_nl_points(C1)
    gm_C2::Int    = add_nl_points(C2)
    gm_P1::Int    = add_nl_points(P1)
    gm_P2::Int    = add_nl_points(P2)
    gm_P3::Int    = add_nl_points(P3)
    # Making Unit Tip Curve
    arc_11   = gmsh.model.occ.add_circle_arc(gm_P1, gm_C1, gm_P2)
    arc_21   = gmsh.model.occ.add_circle_arc(gm_P3, gm_C2, gm_P2)
    arc_12   = gmsh.model.occ.copy([(1, arc_11),])
    arc_22   = gmsh.model.occ.copy([(1, arc_21),])
    gmsh.model.occ.mirror(arc_12, 0.0, 1.0, 0.0, -C[2])
    gmsh.model.occ.mirror(arc_22, 0.0, 1.0, 0.0, -C[2])
    ref_tip_arc_dim_tags = [arc_22[1], arc_12[1], (1, arc_11), (1, arc_21),]
    gmsh.model.occ.rotate(ref_tip_arc_dim_tags, C[1], C[2], C[3], normal_vec[1], normal_vec[2], normal_vec[3], theta)
    #
    n_tip_tags::Vector{Int32} = Int32[i[2] for i in ref_tip_arc_dim_tags]
    for i in 2:num_lobes
        θ_i = 2.0 * α * (i - 1)
        ith_tip_arcs = gmsh.model.occ.copy(ref_tip_arc_dim_tags)
        gmsh.model.occ.rotate(ith_tip_arcs, C[1], C[2], C[3], normal_vec[1], normal_vec[2], normal_vec[3], θ_i)
        append!(n_tip_tags, [i[2] for i in ith_tip_arcs])
    end
    nLobeShape_wire_tag = gmsh.model.occ.add_curve_loop(n_tip_tags)
    nLobeShape_disc_tag = gmsh.model.occ.add_plane_surface([nLobeShape_wire_tag,])
    gmsh.model.occ.synchronize()
    boundary_dimTags = gmsh.model.get_boundary([(2, nLobeShape_disc_tag)], true, false)
    all_curve_tags = gmsh.model.get_entities(1)
    duplicate_curves = [i for i in all_curve_tags if !(i in boundary_dimTags)]
    gmsh.model.occ.remove(duplicate_curves, true)
    gmsh.model.occ.synchronize()
    return nLobeShape_disc_tag
end


# ==========================================
#          CLUSTERS OF REGULAR 2D SHAPES
# ==========================================

function add_circular_discs(
    xyr::Matrix{Float64},
    z_min::Float64,
)::Vector{Int}
    return [gmsh.model.occ.add_disk(ax, ay, z_min, ar, ar) for (ax, ay, ar) in eachrow(xyr)]
end

function add_elliptical_discs(
    xyt_ab::Matrix{Float64},
    z_min::Float64,
)::Vector{Int}
    @assert size(xyt_ab)[2] == 5 "fibre data matrix should contain five columns with xy coordinates of ellipse centre, theta, semi major and minor legths respectively"
    return [make_ellipse((axc, ayc, z_min), ath, aa, ab) for (axc, ayc, ath, aa, ab) in eachrow(xyt_ab)]
end


function add_capsular_discs(
    xyt_ab::Matrix{Float64},
    z_min::Float64,
)::Vector{Int}
    @assert size(xyt_ab)[2] == 5 "fibre data matrix should contain five columns with xy coordinates of capsule centre, theta, semi major and minor legths respectively"
    return [make_capsule((axc, ayc, z_min), ath, aa, ab) for (axc, ayc, ath, aa, ab) in eachrow(xyt_ab)]
end


function add_rectangular_discs(
    xyt_abr::Matrix{Float64},
    z_min::Float64,
)::Vector{Int}
    @assert size(xyt_abr)[2] == 6 "fibre data matrix should contain six columns with xy coordinates of rectangle centre, theta, semi major and minor legths, corner radius respectively"
    return [make_rectangle((ax, ay, z_min), ath, aa, ab, ar) for (ax, ay, ath, aa, ab, ar) in eachrow(xyt_abr)]
end


function add_rpolygons(
    xyt_arn::Matrix{Float64},
    z_min::Float64,
)::Vector{Int}
    @assert size(xyt_arn)[2] == 6 "fibre data matrix should contain six columns with xy coordinates of rpolygon centre, theta, side legnth, corner radius and number of sides respectively"
    return [make_rpolygon((ax, ay, z_min), ath, aa, ar, an) for (ax, ay, ath, aa, ar, an) in eachrow(xyt_arn)]
end


function add_cshapes(
    xyt_ro_ri_alpha::Matrix{Float64},
    z_min::Float64,
)::Vector{Int}
    @assert size(xyt_ro_ri_alpha, 2) == 6 "For C-shaped fibre, data matrix should contain six columns with xy coordinates of its centre, theta, outer radius, inner radius and included angle respectively"
    return [make_cshape((ax, ay, z_min), ath, aro, ari, a_alpha) for (ax, ay, ath, aro, ari, a_alpha) in eachrow(xyt_ro_ri_alpha)]
end


function add_nlobeshapes(
    xyt_ro_rl_n::Matrix{Float64},
    z_min::Float64,
)::Vector{Int}
    @assert size(xyt_ro_rl_n, 2) == 6 "For n-Lobe shaped fibre, data matrix should contain six columns with xy coordinates of its centre, theta, outer radius, lobe radius and number of lobes respectively"
    return [make_nlobeshape((ax, ay, z_min), ath, aro, arl, convert(Int, an)) for (ax, ay, ath, aro, arl, an) in eachrow(xyt_ro_rl_n)]
end


# ==========================================
#          REGULAR 3D SHAPES
# ==========================================


function make_ellipsoid(
    x::Float64, y::Float64, z::Float64,
    azm::Float64, plr::Float64,
    a::Float64, b::Float64, c::Float64,
)::Int
    dim_tag = (3, gmsh.model.occ.add_sphere(x, y, z, 1.0))
    gmsh.model.occ.dilate(dim_tag, x, y, z, a, b, c)  #CHECK_IT major axis or semi-major axis lengths
    gmsh.model.occ.rotate(dim_tag, x, y, z, 0.0, 0.0, 1.0, azm)
    gmsh.model.occ.rotate(dim_tag, x, y, z, sin(azm), -cos(azm), 0.0, (π * 0.5) - plr)  # CHECK_IT if the orientations are as expected?
    return dim_tag[2]
end



# ==========================================
#          CLUSTERS OF REGULAR 3D SHAPES
# ==========================================


function add_spheres(
    xyzr::Matrix{Float64},
)::Vector{Int}
    return [
        gmsh.model.occ.add_sphere(ax, ay, az, ar) for (
            ax, ay, az, ar) in eachrow(xyzr)
    ]
end


"""
    add_ellipsoids() -> Vector{Int}

ellipsoids data must be `Matrix{Float64}` with
x, y, z, azimuthal angle, polar angle, smjrx, intrx, smnrx axis lengths,
as columns, in order.

"""
function add_ellipsoids(
    xyz_ap_abc::Matrix{Float64}
)
    @assert size(xyt_ab, 2) == 8 "Ellipsoids data matrix should contain 8 columns with,
     xyz coordinates of ellipsoid centre,
     azimuthal angle, polar angle, semi major, intermediate and semi minor legths respectively
     BUT $(size(xyt_ab, 2)) is found."
    return [make_ellipsoid(axc, ayc, azc, aazm, aplr, aa, ab, ac) for (
        axc, ayc, azc, aazm, aplr, aa, ab, ac
    ) in eachrow(xyz_ap_abc)]
end






