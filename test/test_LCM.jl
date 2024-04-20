"""
    Testing the Lattice Cell Modelling
"""

using UnitCellModelling

working_folder = mkpath(joinpath(homedir(), "__scrap"))
println("\t\t Testing Lattice Cell Modelling...! \n")

lc_bounds = (-10.0, -10.0, -10.0, 10.0, 10.0, 10.0)
link_eq_r = 2.0
lc = UnitCellModelling.LatticeCell(
    UnitCellModelling.BBox3D(lc_bounds...),
    link_eq_r,
    "Simple_Cubic",
    [0.0, 0.0, 0.0],
)

lc_model_data = make_lattice_cell_model(
    lc;
    small_parameter=1e-06,
    geom_export_paths=(),
    mesh_periodicity=true,
    element_types=("C3D4",),
    min_ele_size_factor=0.01,
    max_ele_size_factor=0.02,
    mesh_opt_algorithm="Netgen",
    show_model=false
)

