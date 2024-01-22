import numpy as np

class Material:
    def __init__(self, **kwargs):
        self.dim = kwargs.get("dim",2)
        self.elastic = kwargs.get("elastic",True)
        self.isotropic = kwargs.get("isotropic",True)
        self.anisotropic = kwargs.get("anisotropic",False)
        self.C_const = kwargs.get("C_const")
        if self.dim == 2:
            self.PlaneStress = kwargs.get("PlaStr",True)
            self.PlaneStrain = kwargs.get("PlaStra",False)
        if self.C_const:
            self.C_const = kwargs.get("C_const")
            self.elastic = False
            self.anisotropic = True
        self.E = kwargs.get("E")
        self.nu = kwargs.get("nu")
        
        self.D = self.Stiffness_Matrix()
    
    def Stiffness_Matrix(self):
        if self.elastic:
            E  = self.E
            nu = self.nu
            if self.dim == 2:
                if self.PlaneStress:
                    D = np.array([[1,  nu,       0],
                                  [nu,  1,       0],
                                  [ 0,  0,(1-nu)/2]])
                    D = E/(1-nu**2)*D
                elif self.PlaneStrain:
                    D = np.array([[1-nu,  nu,        0],
                                  [  nu,1-nu,        0],
                                  [   0,  0,(1-2*nu)/2]])
                    D = E/(1+nu)/(1-2*nu)*D
            elif self.dim == 3:
                D = np.array([[1-nu,  nu,  nu,         0,         0,         0],
                              [  nu,1-nu,  nu,         0,         0,         0],
                              [  nu,  nu,1-nu,         0,         0,         0],
                              [   0,   0,   0,(1-2*nu)/2,         0,         0],
                              [   0,   0,   0,         0,(1-2*nu)/2,         0],
                              [   0,   0,   0,         0,         0,(1-2*nu)/2]])
                D = E/(1+nu)/(1-2*nu)*D
        elif self.anisotropic:
            if len(self.C_const) == 4:
                C11 = self.C_const[0]
                C12 = self.C_const[1]
                C22 = self.C_const[2]
                C66 = self.C_const[3]
                D = np.array([[C11,C12,  0],
                              [C12,C22,  0],
                              [0  ,0  ,C66]])
            elif len(self.C_const) == 9:
                C11 = self.C_const[0]
                C12 = self.C_const[1]
                C13 = self.C_const[2]
                C22 = self.C_const[3]
                C23 = self.C_const[4]
                C33 = self.C_const[5]
                C44 = self.C_const[6]
                C55 = self.C_const[7]
                C66 = self.C_const[8]
                D = np.array([[C11,C12,C13,  0,  0,  0],
                              [C12,C22,C23,  0,  0,  0],
                              [C13,C23,C33,  0,  0,  0],
                              [0  ,0  ,0  ,C44,  0,  0],
                              [0  ,0  ,0  ,  0,C55,  0],
                              [0  ,0  ,0  ,  0,  0,C66]])

        return D

if __name__ == "__main__":
    np.set_printoptions(precision=2, suppress=True)
    # Elastic 2D
    Mat = Material(E=1e3,nu=0.3)
    print('{0:-^60}'.format('Elastic 3D') + '\n' + 
          str(Mat.D))
    # Elastic 3D
    Mat = Material(E=1e3,nu=0.3,dim=3)
    print('{0:-^60}'.format('Elastic 3D') + '\n' + 
          str(Mat.D))
    # Anisotropic 
    C11 = 1e6
    C12 = 0.3e6
    C13 = 0.3e6
    C22 = 1e6
    C23 = 0.3e6
    C33 = 1e6
    C44 = 0.3e6
    C55 = 0.3e6
    C66 = 1e6
    # 2D
    Mat = Material(C_const= [C11,C12,C22,C66],dim=2)
    print('{0:-^60}'.format('Anisotropic 2D') + '\n' + 
          str(Mat.D))
    # 3D
    Mat = Material(C_const= [C11,C12,C13,C22,C23,C33,C44,C55,C66],dim=3)
    print('{0:-^60}'.format('Anisotropic 3D') + '\n' + 
          str(Mat.D))
    