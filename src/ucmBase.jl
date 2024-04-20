udc_inclusions_library = Set(("CIRCLE", "CAPSULE", "RECTANGLE", "ELLIPSE", "REGULARPOLYGON", "CSHAPE", "NLOBESHAPE"))
prc_inclusions_library = Set(("SPHERE", "SPHERO_CYLINDER", "BOX", "ELLIPSOID", "SPHEROID", "OBLATE_SPHEROID", "PROLATE_SPHEROID"))


@with_kw struct BBox2D <: FieldVector{4,Float64}
    xlb::Float64 = -1.0
    ylb::Float64 = -1.0
    xub::Float64 = 1.0
    yub::Float64 = 1.0
    @assert xub >= xlb "While creating 2D bounding box:
     x_max(=$xub) < x_min(=$xlb) is found, please ensure x_max > x_min."
    @assert yub >= ylb "While creating 2D bounding box:
     y_max < y_min is found, please ensure y_max > y_min."
end


@with_kw struct BBox3D <: FieldVector{6,Float64}
    xlb::Float64 = -1.0
    ylb::Float64 = -1.0
    zlb::Float64 = -1.0
    xub::Float64 = 1.0
    yub::Float64 = 1.0
    zub::Float64 = 1.0
    @assert xub >= xlb "While creating 3D bounding box:
     x_max < x_min is found, please ensure x_max > x_min."
    @assert yub >= ylb "While creating 3D bounding box:
     y_max < y_min is found, please ensure y_max > y_min."
    @assert zub >= zlb "While creating 3D bounding box:
     z_max < z_min is found, please ensure z_max > z_min."
end


buffer_bbox(bb::BBox2D, buffer::Float64) = bb .+ (buffer .* (-1.0, -1.0, 1.0, 1.0))
area(bbox::BBox2D) = (bbox.xub - bbox.xlb) * (bbox.yub - bbox.ylb)
buffer_bbox(bb::BBox3D, buffer::Float64) = bb .+ (buffer .* (-1.0, -1.0, -1.0, 1.0, 1.0, 1.0))
volume(bbox::BBox3D) = (bbox.xub - bbox.xlb) * (bbox.yub - bbox.ylb) * (bbox.zub - bbox.zlb)

side_lengths(bbx::BBox2D) = (bbx.xub - bbx.xlb, bbx.yub - bbx.ylb)
side_lengths(bbx::BBox3D) = (bbx.xub - bbx.xlb, bbx.yub - bbx.ylb, bbx.zub - bbx.zlb)


abstract type AbstractUnitCell end

@with_kw struct UDC2D <: AbstractUnitCell
    bbox::BBox2D
    inclusions::Dict{String, Matrix{Float64}}
    @assert issubset(uppercase.(keys(inclusions)), udc_inclusions_library)
end


@with_kw struct UDC3D <: AbstractUnitCell
    bbox::BBox3D
    inclusions::Dict{String, Matrix{Float64}}
    @assert issubset(uppercase.(keys(inclusions)), udc_inclusions_library) "Given inclusion shapes $(uppercase.(keys(inclusions))) doesn't belong to  desired $(udc_inclusions_library)"
end


@with_kw struct PRC <: AbstractUnitCell
    bbox::BBox3D
    inclusions::Dict{String, Matrix{Float64}}
    @assert issubset(uppercase.(keys(inclusions)), prc_inclusions_library)
end

abstract type AbstractLatticeCell <: AbstractUnitCell end

# @with_kw struct LatticeCell <: AbstractUnitCell
#     bbox::BBox3D
#     link_eq_rad::Float64
#     id::String
#     centre::Vector{Float64}
#     @assert link_eq_rad > 0 "Link equivalent radius must be positive."
# end

@with_kw struct SimpleCubicLattice <: AbstractLatticeCell
    bbox::BBox3D
    link_eqr_x::Float64
    link_eqr_y::Float64
    link_eqr_z::Float64
    id::String
    centre::Vector{Float64}
    @assert link_eqr_x > 0 && link_eqr_y > 0 && link_eqr_z > 0 "Link equivalent radius must be positive."
end

const UDC = Union{UDC2D, UDC3D}
const UnitCell3D = Union{UDC3D, PRC, AbstractLatticeCell}
const UNIT_CELL = Union{UDC, UnitCell3D}


side_lengths(uc::AbstractUnitCell) = side_lengths(uc.bbox)

dimension(uc::AbstractUnitCell) = begin
    if isa(uc, UDC2D)
        return 2
    elseif isa(uc, UnitCell3D)
        return 3
    elseif isa(uc, AbstractLatticeCell)
        return 3
    end
end


function buffer_bbox(
    uc::UnitCell3D,
    identifier::Symbol=:ALL,
    ϵ::Float64=1e-06, 
)
    xlb, ylb, zlb, xub, yub, zub = uc.bbox
    if identifier == :ALL
        return BBox3D(xlb-ϵ, ylb-ϵ, zlb-ϵ, xub+ϵ, yub+ϵ, zub+ϵ)
    elseif identifier == :XLB
        return BBox3D(xlb-ϵ, ylb-ϵ, zlb-ϵ, xlb+ϵ, yub+ϵ, zub+ϵ)
    elseif identifier == :YLB
        return BBox3D(xlb-ϵ, ylb-ϵ, zlb-ϵ, xub+ϵ, ylb+ϵ, zub+ϵ)
    elseif identifier == :ZLB
        return BBox3D(xlb-ϵ, ylb-ϵ, zlb-ϵ, xub+ϵ, yub+ϵ, zlb+ϵ)
    elseif identifier == :XUB
        return BBox3D(xub-ϵ, ylb-ϵ, zlb-ϵ, xub+ϵ, yub+ϵ, zub+ϵ)
    elseif identifier == :YUB
        return BBox3D(xlb-ϵ, yub-ϵ, zlb-ϵ, xub+ϵ, yub+ϵ, zub+ϵ)
    elseif identifier == :ZUB
        return BBox3D(xlb-ϵ, ylb-ϵ, zub-ϵ, xub+ϵ, yub+ϵ, zub+ϵ)
    # aϵ = abs(ϵ)
    # return BBox3D(ucbb.xlb-aϵ, ucbb.ylb-ϵ, ucbb.zlb-ϵ, ucbb.xlb+aϵ, ucbb.yub+ϵ, ucbb.zub+ϵ)
    # return BBox3D(ucbb.xlb-ϵ, ucbb.ylb-aϵ, ucbb.zlb-ϵ, ucbb.xub+ϵ, ucbb.ylb+aϵ, ucbb.zub+ϵ)
    # return BBox3D(ucbb.xlb-ϵ, ucbb.ylb-ϵ, ucbb.zlb-aϵ, ucbb.xub+ϵ, ucbb.yub+ϵ, ucbb.zlb+aϵ)
    # return BBox3D(ucbb.xub-aϵ, ucbb.ylb-ϵ, ucbb.zlb-ϵ, ucbb.xub+aϵ, ucbb.yub+ϵ, ucbb.zub+ϵ)
    # return BBox3D(ucbb.xlb-ϵ, ucbb.yub-aϵ, ucbb.zlb-ϵ, ucbb.xub+ϵ, ucbb.yub+aϵ, ucbb.zub+ϵ)
    # return BBox3D(ucbb.xlb-ϵ, ucbb.ylb-ϵ, ucbb.zub-aϵ, ucbb.xub+ϵ, ucbb.yub+ϵ, ucbb.zub+aϵ)
    # IDEA if required, write bboxes for edges too...!    
    end
end


function buffer_bbox(
    uc::UDC2D,
    identifier::Symbol=:ALL,
    ϵ::Float64=1e-06, 
)
    ucbb = uc.bbox
    aϵ = abs(ϵ)
    if identifier == :ALL
        return uc.bbox .+ [-ϵ, -ϵ, ϵ, ϵ]
    elseif identifier == :XLB
        return BBox2D(ucbb.xlb-aϵ, ucbb.ylb-ϵ, ucbb.xlb+aϵ, ucbb.yub+ϵ)
    elseif identifier == :YLB
        return BBox2D(ucbb.xlb-ϵ, ucbb.ylb-aϵ, ucbb.xub+ϵ, ucbb.ylb+aϵ)
    elseif identifier == :XUB
        return BBox2D(ucbb.xub-aϵ, ucbb.ylb-ϵ, ucbb.xub+aϵ, ucbb.yub+ϵ)
    elseif identifier == :YUB
        return BBox2D(ucbb.xlb-ϵ, ucbb.yub-aϵ, ucbb.xub+ϵ, ucbb.yub+aϵ)
    else
        @warn "Invalid identifier is found in buffer_bbox evaluation"
    end
end



"""
Returns SMatrix of adjacent local node labels for each local node of the
given element type.

Each column holds a local node neighbours in order

"""
function adjacent_local_node_indices(
    el_type::Int
)::SMatrix
    if el_type == 1
        return SMatrix{1,2,Int64}([2 1])
    elseif el_type == 2
        return SMatrix{2,3,Int64}([[2, 3] [3, 1] [1, 2]])
    elseif el_type == 3
        return SMatrix{2,4,Int64}([[2, 4] [3, 1] [4, 2] [1, 3]])
    elseif el_type == 16
        return SMatrix{2,8,Int64}(
            [[5, 8] [5, 6] [6, 7] [7, 8] [1, 2] [2, 3] [3, 4] [4, 1]]
        )
    elseif el_type == 4
        return SMatrix{3,4,Int64}([[2, 3, 4] [3, 1, 4] [1, 2, 4] [1, 2, 3]])
    elseif el_type == 6
        return SMatrix{3,6,Int64}(
            [[2, 3, 4] [1, 3, 5] [1, 2, 6] [1, 5, 6] [2, 4, 6] [3, 4, 5]]
        )
    elseif el_type == 5
        return SMatrix{3,8,Int64}(
            [[2, 4, 5] [3, 1, 6] [7, 2, 4] [1, 8, 3] [1, 8, 6] [2, 5, 7] [3, 6, 8] [4, 5, 7]]
        )
    end
end
