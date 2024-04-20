using UnitCellGenerator
using UnitCellModelling
const UCG = UnitCellGenerator
const UCM = UnitCellModelling

working_folder = mkpath(joinpath(homedir(), "WorkOuts", "scratch"))
println("\t\t Testing Unit Cell Modelling...! \n")

# ============================================
#           UNIT CELL GENERATION
# ============================================
r_mean::Float64 = 2.5
r_std::Float64 = 0.0
rve_size::Float64 = 25.0
#
volume_fractions = Dict{String, Float64}(
    "Circles"       => 0.0,
    "Ellipses"      => 0.0,
    "nLobes"        => 0.5,
    "Reg_Polygons"  => 0.0,
    "Rectangles"    => 0.0,
    "Capsules"      => 0.0,
)

#
ruc_bounds::NTuple{4, Float64} = ((-1.0, -1.0, 1.0, 1.0) .* (0.5 * rve_size))  #TODO getting RUC size
const ruc_info = RUC_data(
    bbox=UCG.BBox2D(ruc_bounds...,),
    ssd_ratio=0.07,
    inclusion_distr=RANDOM,
    periodicity=true,
)
# ===================
#   Inclusions data
# ===================
const circ_inc = Inclusion_data(
    volume_fraction=volume_fractions["Circles"],
    shape=Circle,
    size_params=Dict(:RADIUS => Normal(r_mean, r_std),),
)
const caps_inc = Inclusion_data(
    volume_fraction=volume_fractions["Capsules"],
    shape=Capsule,
    size_params=Dict(:SMJRX => Normal(1.42, 0.0), :SMNRX => Normal(0.925, 0.0),),
)
const lobular_inc = Inclusion_data(
    volume_fraction=volume_fractions["nLobes"],
    shape=nLobeShape,
    size_params=Dict(:NLOBES => 4, :EQRAD => Normal(r_mean, r_std), :LOBE_DIST_FACTOR => 0.1,)
)

const ell_inc = Inclusion_data(
    volume_fraction=volume_fractions["Ellipses"],
    shape=Ellipse,
    size_params=Dict(:SMJRX => Normal(2.0, 0.0), :SMNRX => Normal(1.0, 0.0),),
)
const rect_inc = Inclusion_data(
    volume_fraction=volume_fractions["Rectangles"],
    shape=Rectangle,
    size_params=Dict(:SMJRX => Normal(2.0, 0.0), :SMNRX => Normal(1.0, 0.0), :CRAD => Normal(0.2, 0.0),),
)
inclusions_data = generate_unit_cell(ruc_info, (lobular_inc, ell_inc, rect_inc), adjust_ruc_bbox=true,)


# ============================================
#       MODELLING and MESHING in GMESH
# ============================================
rve_bounds =(
    inclusions_data["bbox"][1], inclusions_data["bbox"][2], -0.5*r_mean,
    inclusions_data["bbox"][3], inclusions_data["bbox"][4], 0.5*r_mean,
)
rve_inclusions_data = Dict(k => v for (k,v) in inclusions_data if k!="bbox")
udc_3d = UCM.UDC3D(UCG.BBox3D(rve_bounds...,), rve_inclusions_data,)

ruc_model_data = make_unit_cell_model(
    udc_3d,
    mesh_periodicity=true,
    element_types=(:C3D6, :C3D8),
    geom_export_paths=(),
    extr_dir_num_ele=Int64[3,],
    extr_dir_cum_heights=Float64[1.0,],
    extr_dir_recombine_ele=true,
    min_ele_size_factor=1/15,  #FIXME
    max_ele_size_factor=1/10,   
    mesh_opt_algorithm="Netgen",
    show_mesh_stats=true,
    show_rve=true,
    node_renum_algorithm="",
)
#
# fields of model can be accessed as follows
all_ntags = ruc_model_data["mesh_data"]["all_node_tags"]
all_ncoor = ruc_model_data["mesh_data"]["all_node_coordinates"]
matrix_ele_conn = ruc_model_data["mesh_data"]["matrix_element_connectivity"]
inclusions_ele_conn = ruc_model_data["mesh_data"]["inclusions_element_connectivity"]
total_num_nodes = ruc_model_data["mesh_data"]["total_nodes"]
total_num_3D_elements = ruc_model_data["mesh_data"]["num_3D_elements"]
total_num_2D_elements = ruc_model_data["mesh_data"]["num_2D_elements"]

#
# all_ele_connc = UnitCellModelling.merge_ele_conn([matrix_ele_conn, inclusions_ele_conn])
# println(
#     "Bandwidth before RCM renumbering: ",
#     UnitCellModelling.bandwidth(all_ele_connc, dof=1)
# )
# new_nt_order = UnitCellModelling.update_node_labels([matrix_ele_conn, inclusions_ele_conn])
# println("Updating elem conn ....!")
# all_ele_connc = UnitCellModelling.update_element_connectivity(all_ele_connc, new_nt_order)
# println(
#     "Bandwidth after RCM renumbering: ",
#     UnitCellModelling.bandwidth(all_ele_connc, dof=1)
# )