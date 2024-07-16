from src.Analysis import model
from src.Material import material
from src.Assemble import assemble
from src.Boundary import Boundary
from src.Solver import solver
from src.PostProcess import Plot2D
import matplotlib.pyplot as plt
import meshio
import numpy as np

# *****************************************************
#  Create an instance of the model
# *****************************************************
model1 = model(NDof=2)
# *****************************************************
#  Create mesh (exported from Abaqus)
# *****************************************************
# Read the mesh information
mesh = meshio.read('Job-1.inp')
# Read the node set information
nset_top = mesh.point_sets['Top']
nset_bottom = mesh.point_sets['Bottom']
nset_left = mesh.point_sets['Left']
nset_right = mesh.point_sets['Right']
# Add the mesh to the model
model1.mesh_add(mesh)
# *****************************************************
#  Create the material
# *****************************************************
# Define the material properties
mat = material([210e9, 0.3], name='steel', type='isotropic')
# *****************************************************
#  Create boundary conditions
# *****************************************************
boundary = Boundary()
boundary.set_bc(nset_bottom, 1, 0)
boundary.set_bc(nset_bottom, 2, 0)
boundary.set_bc(nset_top, 2, 0.05)
# *****************************************************
#  Assemble the global stiffness matrix
# *****************************************************
assemble1 = assemble()
# *****************************************************
#  Solve the system of equations
# *****************************************************
Result = solver(steps=10)
# Plot2D(mesh, disp.disp,deformed=True,wireframe=True)
disp_vec = Result.disp
disp_vec = disp_vec.astype(np.float32)
ux = disp_vec[0::2]
uy = disp_vec[1::2]
# *****************************************************
#  Post-process
# *****************************************************
Plot2D(mesh, ux, deformed=False, wireframe=True)
# *****************************************************
#  Write the results to a VTK file
# *****************************************************
mesh.point_data['ux'] = ux
mesh.point_data['uy'] = uy
mesh.point_data['u'] = np.linalg.norm([ux, uy], axis=0)
meshio.write_points_cells('Job-1.vtk', mesh.points, mesh.cells, point_data=mesh.point_data)