import numpy as np
from src.Analysis import model



class Boundary:
    Penalty = 1e16
    Force = None
    def __init__(self):
        
        pass
    
    def set_bc(self, node_set, degree, value):
        """Set the boundary condition
        
        Args:
            node_set (list): list of node numbers
            degree (str): degree of freedom
            value (float): value of the boundary condition
        """
        # Update boundary information
        mesh = model.mesh
        NDof = model.NDof
        # Find node set name
        for key, point_set in mesh.point_sets.items():
            if np.array_equal(node_set, point_set):
                set_name = key+'_'+str(degree)
                break
        # Add boundary condition to the model
        model.bc_cell_set_add(set_name, node_set, degree, value)
        # Calculate the force vector
        Force = np.zeros((len(mesh.points)*NDof))
        Force[NDof*node_set+degree-1] = value * Boundary.Penalty
        # Add the force to the global force vector
        if Boundary.Force is None:
            Boundary.Force = Force
        Boundary.Force += Force
    