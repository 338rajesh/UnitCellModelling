# Unit Cell Model

This module is intended to perform geometric modelling and meshing of an Unit Cell using [`gmsh`](http://gmsh.info/), an open source code.

It exports a single function `create()` which takes the following two arguments.

 1. `unit_cell::AbstractUnitCell`.

    * `UDC2D <: AbstractUnitCell`
    * `UDC3D <: AbstractUnitCell`
    * `PRC <: AbstractUnitCell`
    * `UDC = Union{UDC2D, UDC3D}`
    * `UnitCell3D = Union{UDC3D, PRC}`

 2. `options::Union{Tuple,Dict{Symbol,Any}} = ()`, By default the following key-value pairs will be used as options

    ```julia
    :small_parameter => 1e-06,
    :geom_export_paths => (joinpath(homedir(), "unit_cell.step"),),
    :extrusion_direction => "XY->Z",
    :mesh_periodicity => true,
    :element_types => (:DEFAULT,),
    :num_ele_in_extr_dir => Int64[],
    :cum_heights => Float64[],
    :recombine_ele_in_extr_dir => true,
    :mesh_min_size_factor => 2.0,
    :mesh_max_size_factor => 1.0,
    :mesh_opt_algorithm => "Netgen",
    :show_mesh_stats => true,
    :show_rve => true,
    ```

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
