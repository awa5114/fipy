from fipy import CylindricalGrid2D, Grid2D
import pyvista


mesh = Grid2D(dx = 1.0, dy = 1.0, nx = 2, ny = 2)

ugrid= pyvista.UnstructuredGrid(mesh.VTKCellDataSet._vtk_obj)
ugrid.save('broken_mesh.vtk', binary=False)


plotter = pyvista.Plotter()
plotter.set_background('white')
plotter.add_mesh(ugrid, style='wireframe', color='black')
plotter.add_bounding_box(color='red')

plotter.view_xy()
plotter.show()

