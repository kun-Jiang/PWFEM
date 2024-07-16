import numpy as np
from src.Analysis import model
from src.Boundary import Boundary

def gauss_quad_2d(n):
    if n == 4:
        gauss_points = np.array([
            [-1/np.sqrt(3), -1/np.sqrt(3)],
            [1/np.sqrt(3), -1/np.sqrt(3)],
            [1/np.sqrt(3), 1/np.sqrt(3)],
            [-1/np.sqrt(3), 1/np.sqrt(3)]
        ])
        gauss_weights = np.array([1, 1, 1, 1])
        return gauss_points, gauss_weights
    
def shape_func(xi, eta):
    # Shape functions for a 4-node quadrilateral element
    N1 = 0.25 * (1 - xi) * (1 - eta)
    N2 = 0.25 * (1 + xi) * (1 - eta)
    N3 = 0.25 * (1 + xi) * (1 + eta)
    N4 = 0.25 * (1 - xi) * (1 + eta)
    # The derivative of the shape functions w.r.t. xi and eta
    dN_dx = np.array([
        [-0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta)],
        [-0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi)]
    ])
    return np.array([N1, N2, N3, N4]), dN_dx

class assemble:
    K = None
    def __init__(self):
        self.mesh = model.mesh
        self.NDof = model.NDof
        self.NNodes = model.NNodes
        self.nodes = self.mesh.points
        self.elements = self.mesh.cells_dict['quad']
        D = model.mat_cell_sets['steel']['matrix']
        assemble.K = self.assemble_stiffness(D)
        
    def assemble_stiffness(self,D):
        # Loop through the elements
        K = np.zeros((self.NNodes*self.NDof,self.NNodes*self.NDof))
        for element in self.elements:
            connectivity = element.astype(int)
            # Loop through the nodes in the element
            # coords = np.zeros((4, 2))
            coords = self.nodes[connectivity]
            gauss_points,gauss_weight = gauss_quad_2d(4)
            for m in range(len(gauss_points)):
                xi, eta = gauss_points[m]
                weight = gauss_weight[m]
                # Calculate the shape function
                N,dN_dx = shape_func(xi,eta)
                J = np.zeros((2, 2))
                for k in range(4):
                    J[0, 0] += dN_dx[0, k] * coords[k, 0]
                    J[0, 1] += dN_dx[0, k] * coords[k, 1]
                    J[1, 0] += dN_dx[1, k] * coords[k, 0]
                    J[1, 1] += dN_dx[1, k] * coords[k, 1]
                detJ = np.linalg.det(J)
                # Assembly the stiffness matrix
                for i,Int_i in enumerate(connectivity):
                    B_i = np.array([
                        [dN_dx[0,i],       0.0],
                        [0.0,       dN_dx[1,i]],
                        [dN_dx[1,i],dN_dx[0,i]]
                    ])
                    for j,Int_j in enumerate(connectivity):
                        B_j =  np.array([
                            [dN_dx[0,j],       0.0],
                            [0.0,       dN_dx[1,j]],
                            [dN_dx[1,j],dN_dx[0,j]]
                        ])
                        C = np.dot(np.dot(B_j.T,D),B_i)
                        K[2*Int_i:2*Int_i+2,2*Int_j:2*Int_j+2] += C*detJ*weight
        # Add penalty terms
        bc_cell_sets = model.bc_cell_set_get()
        for key, bc_set in bc_cell_sets.items():
            node_set = bc_set[key]
            dof = bc_set['dof']
            K[self.NDof*node_set+dof-1,self.NDof*node_set+dof-1] += Boundary.Penalty
        return K