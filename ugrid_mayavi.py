from numpy import array, random
from tvtk.api import tvtk
from mayavi import mlab

points = array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],  # tets
                [1, 0, 0], [2, 0, 0], [1, 1, 0], [1, 0, 1],
                [2, 0, 0], [3, 0, 0], [2, 1, 0], [2, 0, 1],
                ], 'f')
tets = array([[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]])
tet_type = tvtk.Tetra().cell_type
ug = tvtk.UnstructuredGrid(points=points)
ug.set_cells(tet_type, tets)



fig = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0),
                  figure=ug.class_name[3:])
surf = mlab.pipeline.surface(ug, opacity=0.1)
mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                        color=(0, 0, 0), )
mlab.show()