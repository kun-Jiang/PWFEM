import numpy as np
from src.Analysis import model




class material():
    def __init__(self, property=None, name=None, type=None):
        """Define a material class

        Args:
            paras (list): material properties
            name (str): material name
            type (_type_): mateiral type: isotropic, orthotropic, anisotropic
        """
        self.property = property
        self.name     = name
        self.type     = type
        self.name     = name
        self.D        = None
        self.add_material()
        
    def add_material(self):
        mat_cell_sets = model.mat_cell_set_get()
        if self.property is not None:
            self.D = self.cal_stiffness_matrix()
        if self.name is None:
            self.name = "Material "+str(len(mat_cell_sets)+1)
            if self.name in mat_cell_sets.keys():
                self.name = self.name + "_1"

            
        
        model.mat_cell_set_add(cell_name=self.name,
                           cell_set={'mat':self.property,'matrix':self.D,'type':self.type})
        
    def cal_stiffness_matrix(self):
        E = self.property[0]
        nu = self.property[1]
        
        D = np.array([[1,  nu,       0],
                      [nu,  1,       0],
                      [ 0,  0,(1-nu)/2]])
        D = E/(1-nu**2)*D
        return D
        