using gmsh
using UnitCellModelling

gmsh.initialize()
gmsh.model.add("test")

test_cshape = false
test_nlobe_shape = true

if test_nlobe_shape
    UnitCellModelling.make_nlobeshape(
        (0.0, 0.0, 0.0),
        π * 0.0, 20.0, 10.0, 3,
    )
elseif test_cshape
    UnitCellModelling.make_cshape(
        (0.0, 0.0, 0.0),
        0.25*π,
        20.0, 15.0, π*0.999
    )
end

gmsh.model.occ.synchronize()
gmsh.fltk.run()
gmsh.finalize()
