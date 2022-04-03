using gmsh
using UnitCellModelling

gmsh.initialize()
gmsh.model.add("test")

test_cshape = false
test_nlobe_shape = true


nls_tag = UnitCellModelling.make_nlobeshape(
    (0.0, 0.0, 0.0),
    π * 0.0, 10.0, 5.0, 3,
)

csh_tag = UnitCellModelling.make_cshape(
    (15.0, 10.0, 0.0),
    1.25*π,
    5.0, 2.0, π*0.999
)

rect = gmsh.model.occ.add_rectangle(-25.0, -25.0, 0.0, 50.0, 50.0)
frags, _ = gmsh.model.occ.fragment([(2, rect)], [(2, csh_tag), (2, nls_tag)])
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
gmsh.option.set_number("Mesh.MeshSizeMin", 0.5)
gmsh.option.set_number("Mesh.MeshSizeMax", 1.8)
gmsh.model.mesh.generate(2)
gmsh.fltk.run()
gmsh.finalize()
