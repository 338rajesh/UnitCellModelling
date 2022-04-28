
using UnitCellModelling
using FileIO

# ucodb_path = raw"C:\Users\admin\WorkOuts\scratch\field_plots\single_inclusion_34020.jld2"
ucodb_path = raw"C:\Users\admin\WorkOuts\scratch\field_plots\single_inclusion.jld2"

results_dir = raw"C:\Users\admin\WorkOuts\scratch\field_plots"

uc_ODB = load(ucodb_path)



println(typeof(uc_ODB))
println(sizeof(uc_ODB))

elastic_load_cases = [i for i in keys(uc_ODB) if startswith(i, "FVG_eps")]
thermal_load_cases = [i for i in keys(uc_ODB) if startswith(i, "FVG_Delta")]
piezo_electric_load_cases = [i for i in keys(uc_ODB) if (startswith(i, "FVG_eps") && endswith(i, "piezo_electric"))]

UnitCellModelling.plot_nodal_fields(
     uc_ODB,
     elastic_load_cases,
     Dict(
        #  "Sigma 11" => 10,
         "Sigma 22" => 11,
        #  "Sigma 33" => 12,
         "Sigma 23" => 13,
        #  "Sigma 31" => 14,
        #  "Sigma 12" => 15,
     ),
     plots_dir=results_dir,
     field_type="S",
     plot_type=".pos"
 )


# UnitCellModelling.plot_nodal_fields(
#     uc_ODB,
#     thermal_load_cases,
#     Dict(
#         "DeltaT 11" => 2,
#         "DeltaT 22" => 3,
#         "DeltaT 33" => 4,
#         "Heat Flux q11" => 5,
#         "Heat Flux q22" => 6,
#         "Heat Flux q33" => 7,
#     ),
#     plots_dir=results_dir,
#     field_type="S" 
# )

#UnitCellModelling.plot_nodal_fields(
#    uc_ODB,
#    piezo_electric_load_cases,
#  Dict(
#        "Sigma 11" => 14,
#       "Sigma 22" => 15,
#        "Sigma 33" => 16,
#        "Sigma 23" => 17,
#        "Sigma 31" => 18,
#        "Sigma 12" => 19,
#        "Electric Displacements T11" => 20,
#        "Electric Displacements T22" => 21,
#        "Electric Displacements T33" => 22,
#    ),
#    plots_dir=results_dir,
#    field_type="S" 
#)


println("Done! results are saved at $results_dir")
