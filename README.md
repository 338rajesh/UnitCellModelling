# Unit Cell Model

This module is intended to perform geometric modelling and meshing of an Unit Cell using [`gmsh`](http://gmsh.info/), an open source code.

## Installation

```julia
julia> ]add https://github.com/338rajesh/UnitCellModelling.jl 
```

## Usage

It exports a single function `create(unit_cell; <kw_args>)` which takes the following two arguments.

### Positional Argument(s)

1. `unit_cell::AbstractUnitCell`, where `AbstractUnitCell` has three subtypes `UDC2D`, `UDC3D`, `PRC`. Also, another two types `UDC = Union{UDC2D, UDC3D}` and `UnitCell3D = Union{UDC3D, PRC}` for convenience.

### Keyword Arguments

* `small_parameter::Float64 = 1e-06`,

   This parameter is used wherever the geometric precision is considered during the modelling process.

* `geom_export_paths::Tuple{Vararg{String}} = ()`
  
  Unit cell can be exported in various formats such as `.step`, `.brep`,...etc.

* `extrusion_direction::String = "XY->Z"`,

   For example, with `"XY->Z"`, Unit cell is first modelled in `XY` plane and then extruded in `Z` direction.

* `extr_dir_num_ele::Vector{Int64} = Int64[]`,

   Number of elements in the direction of extrusion/thickness. This is applicable only when mesh is generated with extrusion.

* `extr_dir_cum_heights::Vector{Float64} = Float64[]`,

   Cumulative height of elements in the extrusion/thickness direction. For example, with `extr_dir_num_ele=3`, cumulative heights [0.2, 0.6, 1.0] imply that element thickness take 20%, 40% and 40% thickness.

* `extr_dir_recombine_ele::Bool = true`,
* `mesh_periodicity::Bool = true`,
* `element_types::Tuple{Vararg{Symbol}} = (:DEFAULT,)`,
* `min_ele_size_factor::Float64 = 1.0`,
* `max_ele_size_factor::Float64 = 2.0`,
* `mesh_opt_algorithm::String = "Netgen"`,
* `show_mesh_stats::Bool = true`,
* `show_rve::Bool = true`,

### Returns

1. Unit Cell Geometry
   1.
2. mesh_data
   1. Element Tags/Labels and corresponding nodal_connectivity
   2. Node Tags/Labels and corresponding coordinates

As of now, fibres with the following shapes can be modelled,

1. Circles `add_circular_discs()`
2. Regular Polygon `make_rpolygon()` or `add_rpolygons()`
3. Capsule `make_capsule()` or `add_capsular_discs()`
4. Ellipse `make_ellipse()` or `add_elliptical_discs()`
5. Rectangle `make_rectangle()` or `add_rectangular_discs()`
