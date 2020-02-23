from fipy import Grid2D
import pyvista

# Create Mesh
mesh = Grid2D(dx = 1.0, dy = 1.0, nx = 2, ny = 2)
ug = mesh.VTKCellDataSet
#ug.set_cells([7,7,7,7], ug.cell_locations_array.to_array(), ug.get_cells()) # overwrite to
ugrid= pyvista.UnstructuredGrid(ug._vtk_obj)

# Plot Mesh
plotter = pyvista.Plotter()
plotter.set_background('white')
plotter.add_mesh(ugrid, style='wireframe', color='black')
plotter.add_bounding_box(color='red')

plotter.view_xy()
plotter.show()

