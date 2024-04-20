module UnitCellModelling
	__precompile__()
	
	using gmsh
	using Parameters
	using StaticArrays
	using SparseArrays


	include("ucmBase.jl")  # custom type definitions and their methods
	include("gshapes.jl")  # gmsh based shape(s) creation
	include("geometric_modelling.jl")  # create the geometries of the whole RVE/UnitCell
	include("lattice_cell_modelling.jl")  # create the geometries of the lattice cells
	include("meshing.jl")  # Applying periodicity constraints, generating mesh and its statistics
	include("node_labelling.jl")  # node re-numbering by Reverse Cuthill-McKee implementation.
	include("visualization.jl")  # local fields plotting
	include("unit_cell_mesh_model.jl")  # top-level implementation of this package

	export make_udc_model
	export make_lattice_cell_model
	export plot_a_field

end  # of UnitCellModelling module
