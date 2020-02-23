from fipy import CellVariable, Gmsh2D
import pygmsh, pyvista

geom = pygmsh.built_in.Geometry()
p1 = geom.add_point((0.,0.,0.), lcar = 1.0)
p2 = geom.add_point((2.,0.,0.), lcar = 1.0)
p3 = geom.add_point((0.,2.,0.), lcar = 1.0)
p4 = geom.add_point((2.,2.,0.), lcar = 1.0)
l6 = geom.add_line(p1, p2)
l7 = geom.add_line(p2, p4)
l8 = geom.add_line(p4, p3)
l9 = geom.add_line(p3, p1)
ll10 = geom.add_line_loop((l6, l7, l8, l9))
s11 = geom.add_surface(ll10)
geom.set_transfinite_lines([l6, l7, l8, l9], 3)
geom.set_transfinite_surface(s11)
geom.set_recombined_surfaces([s11])
geom_string = geom.get_code()
mesh = Gmsh2D(geom_string)

phi = CellVariable(name = "solution variable", mesh = mesh, value = 0.)

ugrid= pyvista.UnstructuredGrid(mesh.VTKCellDataSet._vtk_obj)
ugrid.save('working_mesh.vtk', binary=False)


plotter = pyvista.Plotter()
plotter.set_background('white')
plotter.show_grid(color='red')
plotter.add_mesh(ugrid, style='wireframe', color='black')
plotter.add_bounding_box(color='red')

plotter.view_xy()

plotter.show()

