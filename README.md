# Unit Cell Modelling

Geometric modelling of a unit cell and its meshing can be obtained using this package. This task is accomplished using an open source code, [`gmsh`](http://gmsh.info/) in the background. 

> **Note**: Here, unit cell implies any spatially repeating rectangular/cuboid domain. It may have inclusions of various shapes and their spatial arrangement may follow any distribution.

---

## Installation

```julia
julia> ]add https://github.com/338rajesh/UnitCellModelling.jl 
```

---

## Usage

This module exports the following functions.

```julia
make_unit_cell_model(unit_cell; <kwargs>)
plot_nodal_fields(uc_odb, load_case_IDs, nodal_fvariables; <kwargs>)
```

- `unit_cell::AbstractUnitCell`
  
   Defining an unit cell requires a bbox and inclusion geometry data. This information is supplied to modelling process using this `unit_cell` argument. For example, the following lines define a 3D unit cell of unidirectional composite with circular fibres.

   ```julia
   julia> using UnitCellModelling
   julia> uc_bbox = UnitCellModelling.BBox3D(-1.0,-1.0,-0.2, 1.0, 1.0, 0.2)
   julia> uc_incl_data = Dict("CIRCLE" => [0.0 0.0 0.15; 0.5, 0.5 0.5 0.2])
   julia> unit_cell = UnitCellModelling.UDC3D(uc_bbox, uc_incl_data)
   ```

   Note, from the BBox3D definition, that the unit cell dimensions are `2.0, 2.0, 0.4` in `x, y, z` directions and its centre is at the origin. Also, two circles with radii `0.15` and `0.2` are placed at `(0.0, 0.0)` and `(0.5, 0.5)`

   `AbstractUnitCell` has three subtypes `UDC2D`, `UDC3D`, `PRC`. Also, another two types `UDC = Union{UDC2D, UDC3D}` and `UnitCell3D = Union{UDC3D, PRC}` are defined for convenience.

   [//]: # (write details of possible fibre shapes...etc)

- kwargs, Keyword Arguments
  - `small_parameter::Float64 = 1e-06`,

      this parameter is used wherever the geometric precision is considered during the modelling process.

  - `geom_export_paths::Tuple{Vararg{String}} = ()`,

      Unit cell can be exported in various formats such as `.step`, `.brep`,...etc.

  - `extrusion_direction::String = "XY->Z"`,

     For example, with `"XY->Z"`, Unit cell is first modelled in `XY` plane and then extruded in `Z` direction.

  - `extr_dir_num_ele::Vector{Int64} = Int64[]`,

     Number of elements in the direction of extrusion/thickness. This is applicable only when mesh is generated with extrusion.

  - `extr_dir_cum_heights::Vector{Float64} = Float64[]`,

     Cumulative height of elements in the extrusion/thickness direction. For example, with `extr_dir_num_ele=3`, cumulative heights [0.2, 0.6, 1.0] imply that element thickness take 20%, 40% and 40% thickness.

  - `extr_dir_recombine_ele::Bool = true`,
  - `mesh_periodicity::Bool = true`,
  - `element_types::Tuple{Vararg{Symbol}} = (:DEFAULT,)`,
  - `min_ele_size_factor::Float64 = 1.0`,
  - `max_ele_size_factor::Float64 = 2.0`,
  - `mesh_opt_algorithm::String = "Netgen"`,
  - `show_mesh_stats::Bool = true`,
  - `node_renum_algorithm = "RCM"`,
  - `show_rve::Bool = true`,

- Returns

   A dictionary with the following key-value pairs.

  - `"mesh_data"` being the key, value is a dictionary with the following key-value pairs,
    - `"all_node_tags"` => `Vector{Int}`
    - `"all_node_coordinates"` => `Matrix{Float64}`
    - `"matrix_element_connectivity"` => `Matrix{Int}`
    - `"inclusions_element_connectivity"` => `Matrix{Int}`
    - `"total_nodes"` => Int
    - `"num_2D_elements"` => Int
    - `"num_3D_elements"` => Int

### Acceptable inclusion shapes
+ 2D shapes:
   + `"CIRCLE"`
   + `"CAPSULE"`
   + `"RECTANGLE"`
   + `"ELLIPSE"`
   + `"REGULARPOLYGON"`
   + `"CSHAPE"`
   + `"NLOBESHAPE"`

---

As of now, fibres with the following shapes can be modelled,

1. Circles `add_circular_discs()`
2. Regular Polygon `make_rpolygon()` or `add_rpolygons()`
3. Capsule `make_capsule()` or `add_capsular_discs()`
4. Ellipse `make_ellipse()` or `add_elliptical_discs()`
5. Rectangle `make_rectangle()` or `add_rectangular_discs()`
