import meshio
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import griddata

# This is a simple script to solve linear elastic material based 
# on the finite element method.
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

if __name__ == '__main__':
    # Read the mesh information exported from Abaqus
    mesh = meshio.read('Job-1.inp')
    nodes = mesh.points
    elements = mesh.cells_dict['quad']
    nset_top = mesh.point_sets['Top']
    nset_bottom = mesh.point_sets['Bottom']
    nset_left = mesh.point_sets['Left']
    nset_right = mesh.point_sets['Right']
    # Define the degree of freedom
    dof = ['x','y']
    # Define the material properties
    Eyoung = 1e9
    Poisson = 0.3
    D = np.array([
        [1.0,   Poisson,    0.0],
        [Poisson,   1.0,    0.0],
        [0.0,       0.0,    (1-Poisson)/2.0]
    ])*Eyoung/(1-np.power(Poisson,2))
    # Loop through the elements
    K = np.zeros((len(nodes)*len(dof),len(nodes)*len(dof)))
    for element in elements:
        connectivity = element.astype(int)
        # Loop through the nodes in the element
        # coords = np.zeros((4, 2))
        coords = nodes[connectivity]
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
    # meshio.write('Job-1.vtk', mesh)
    Penalty=1.0e16
    F = np.zeros((len(nodes)*len(dof)))
    F[2*nset_bottom] = 0
    K[2*nset_bottom,2*nset_bottom] += Penalty
    F[2*nset_bottom+1] = 0
    K[2*nset_bottom+1,2*nset_bottom+1] += Penalty
    F[2*nset_top+1] = 0.05*Penalty
    K[2*nset_top+1,2*nset_top+1] += Penalty
    # 增量加载和迭代求解
    disp = np.zeros((len(nodes) * len(dof)))  # 初始位移为零
    load_increment = 10  # 加载步长
    tolerance = 1e-6  # 收敛容限
    max_iter = 50  # 最大迭代次数
    total_load_steps = 10
    for step in range(total_load_steps):
        delta_F = F / total_load_steps  # 每一步的增量力
        residual = delta_F
        iter_count = 0
        while np.linalg.norm(residual) > tolerance and iter_count < max_iter:
            delta_u = np.linalg.solve(K, residual)  # 计算增量位移
            disp += delta_u  # 更新总位移
            # 更新残差
            internal_force = np.dot(K, delta_u)
            residual = (delta_F - internal_force)/Penalty
            iter_count += 1
        if iter_count == max_iter:
            print(f"Warning: Load step {step+1} did not converge")
    print([max(disp),min(disp)])
    fig = plt.figure()
    x,y = nodes[:,0],nodes[:,1]
    X,Y=np.meshgrid(x,y)
    Z=griddata((x,y),disp[0::2],(X,Y),method='linear')
    # plt.contourf(X,Y,Z,vmin=min(disp),vmax=max(disp),levels=20)
    plt.contourf(X,Y,Z,levels=20)
    plt.colorbar()
    for element in elements:
        quad = nodes[element]
        # print(np.append(quad, quad[0:1],axis=0))
        plt.plot(*zip(*np.append(quad, quad[0:1], axis=0)), color='k', linewidth=0.5)
    plt.show()