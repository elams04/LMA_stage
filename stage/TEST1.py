
import gmsh
import sys

gmsh.initialize(sys.argv)
gmsh.model.add("Test3")
lc=0.01
gmsh.model.geo.addPoint(0,0,0,lc,1)
gmsh.model.geo.addPoint(.1, 0, 0, lc, 2)
gmsh.model.geo.addPoint(.1, .3, 0, lc, 3)
p4 = gmsh.model.geo.addPoint(0, .3, 0, 1)


gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(3, 2, 2)
gmsh.model.geo.addLine(3, p4, 3)
gmsh.model.geo.addLine(4, 1, p4)
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()