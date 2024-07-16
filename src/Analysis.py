import numpy as np


class model:
    mat_cell_sets = {}
    bc_cell_sets = {}
    mesh = None
    NDof = int
    NNodes = int
    NElements = int
    def __init__(self,NDof):
        model.NDof = NDof
        # self.material = material()
        # model.cell_sets
        pass
    @classmethod
    def mat_cell_set_get(self):
        return model.mat_cell_sets
    @classmethod
    def mat_cell_set_add(self, cell_name, cell_set):
        model.mat_cell_sets[cell_name] = cell_set
        
    @classmethod
    def bc_cell_set_get(self):
        return model.bc_cell_sets
    @classmethod
    def bc_cell_set_add(self, cell_name, cell_set, degree, value):
        model.bc_cell_sets[cell_name] = {cell_name:cell_set, 
                                         'dof':degree, 
                                         'value':value}
        
    @classmethod
    def mesh_add(self, mesh):
        model.mesh = mesh
        model.NNodes = len(mesh.points)
        model.NElements = len(mesh.cells)