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
    @assert issubset(uppercase.(keys(inclusions)), udc_inclusions_library)
end


@with_kw struct PRC <: AbstractUnitCell
    bbox::BBox3D
    inclusions::Dict{String, Matrix{Float64}}
    @assert issubset(uppercase.(keys(inclusions)), prc_inclusions_library)
end

const UDC = Union{UDC2D, UDC3D}
const UnitCell3D = Union{UDC3D, PRC}

side_lengths(uc::AbstractUnitCell) = side_lengths(uc.bbox)

dimension(uc::AbstractUnitCell) = begin
    if isa(uc, UDC2D)
        return 2
    elseif isa(uc, UnitCell3D)
        return 3
    end
end


function buffer_bbox(
    uc::UnitCell3D,
    identifier::Symbol=:ALL,
    ϵ::Float64=1e-06, 
)
    ucbb = uc.bbox
    aϵ = abs(ϵ)
    if identifier == :ALL
        return uc.bbox .+ [-ϵ, -ϵ, -ϵ, ϵ, ϵ, ϵ]
    elseif identifier == :XLB
        return BBox3D(ucbb.xlb-aϵ, ucbb.ylb-ϵ, ucbb.zlb-ϵ, ucbb.xlb+aϵ, ucbb.yub+ϵ, ucbb.zub+ϵ)
    elseif identifier == :YLB
        return BBox3D(ucbb.xlb-ϵ, ucbb.ylb-aϵ, ucbb.zlb-ϵ, ucbb.xub+ϵ, ucbb.ylb+aϵ, ucbb.zub+ϵ)
    elseif identifier == :ZLB
        return BBox3D(ucbb.xlb-ϵ, ucbb.ylb-ϵ, ucbb.zlb-aϵ, ucbb.xub+ϵ, ucbb.yub+ϵ, ucbb.zlb+aϵ)
    elseif identifier == :XUB
        return BBox3D(ucbb.xub-aϵ, ucbb.ylb-ϵ, ucbb.zlb-ϵ, ucbb.xub+aϵ, ucbb.yub+ϵ, ucbb.zub+ϵ)
    elseif identifier == :YUB
        return BBox3D(ucbb.xlb-ϵ, ucbb.yub-aϵ, ucbb.zlb-ϵ, ucbb.xub+ϵ, ucbb.yub+aϵ, ucbb.zub+ϵ)
    elseif identifier == :ZUB
        return BBox3D(ucbb.xlb-ϵ, ucbb.ylb-ϵ, ucbb.zub-aϵ, ucbb.xub+ϵ, ucbb.yub+ϵ, ucbb.zub+aϵ)
    # IDEA if required, write bboxes for edges too...!    
    end
end

