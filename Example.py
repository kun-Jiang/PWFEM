import pygmsh
import numpy as np
from src.Preprocess import Preprocess
from src.Material import Material as Mat
from src.LinearSolve import LinearSolve as LiSolve
from src.Postprocess import PostProcess


with pygmsh.geo.Geometry() as geom:
    geom.add_polygon(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ],
        mesh_size=2,
        )
    mesh = geom.generate_mesh()
nodes = mesh.points
elements = mesh.cells_dict["triangle"]


# Preprocess.mesh_visual(model)

Mat(E=1e6,nu=0.3)

result = LiSolve(nodes,elements)

s = np.array(result.stress)

PostProcess(nodes,elements,s)

