# Unit Cell Model

This module is intended to perform geometric modelling and meshing of an Unit Cell using [`gmsh`](http://gmsh.info/), an open source code.

It exports a single function `create()` which takes the following two arguments.

## Arguments

1. `unit_cell::AbstractUnitCell`, where
   * `UDC2D <: AbstractUnitCell`
   * `UDC3D <: AbstractUnitCell`
   * `PRC <: AbstractUnitCell`
   * `UDC = Union{UDC2D, UDC3D}`
   * `UnitCell3D = Union{UDC3D, PRC}`

## Keyword Arguments

* `small_parameter::Float64 = 1e-06`,
* `geom_export_paths::Tuple{Vararg{String}} = ()`,
* `extrusion_direction::String = "XY->Z"`,
* `extr_dir_num_ele::Vector{Int64} = Int64[]`,
* `extr_dir_cum_heights::Vector{Float64} = Float64[]`,
* `extr_dir_recombine_ele::Bool = true`,
* `mesh_periodicity::Bool = true`,
* `element_types::Tuple{Vararg{Symbol}} = (:DEFAULT,)`,
* `min_ele_size_factor::Float64 = 1.0`,
* `max_ele_size_factor::Float64 = 2.0`,
* `mesh_opt_algorithm::String = "Netgen"`,
* `show_mesh_stats::Bool = true`,
* `show_rve::Bool = true`,

It is designed to return the following information

1. Unit Cell Geometry
   1. Exporting unit cell as `.step`, `.brep`,...etc. formats.
2. Unit Cell Mesh data
   1. Element Tags/Labels and corresponding nodal_connectivity
   2. Node Tags/Labels and corresponding coordinates

As of now, fibres with the following shapes can be modelled,

1. Circles `add_circular_discs()`
2. Regular Polygon `make_rpolygon()` or `add_rpolygons()`
3. Capsule `make_capsule()` or `add_capsular_discs()`
4. Ellipse `make_ellipse()` or `add_elliptical_discs()`
5. Rectangle `make_rectangle()` or `add_rectangular_discs()`
