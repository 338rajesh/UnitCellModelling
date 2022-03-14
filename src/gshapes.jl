# ==========================================
#          REGULAR 2D SHAPES
# ==========================================


function make_rpolygon(
    xyz_centre::NTuple{3,Float64},
    theta::Float64,
    side_len::Float64,
    rf::Float64,
    num_sides::Real;
    normal_vec::NTuple{3,Float64} = (0.0, 0.0, 1.0)
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
    gmsh.model.occ.rotate(dim_tag, x, y, z, sin(azm), -cos(azm), 0.0, (π*0.5) - plr)  # CHECK_IT if the orientations are as expected?
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






