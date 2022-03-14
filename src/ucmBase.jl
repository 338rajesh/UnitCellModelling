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