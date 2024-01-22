import numpy as np

class LinearSolve:
    def __init__(self,nodes,elements):
        self.nodes = nodes
        self.element = elements
        self.B = self.cal_strain_matix()
        self.d = [0,0,0,-0.1,0.0,0.0]
        self.E = 1e6
        self.nu = 0.3
        self.D = self.cal_D()
        self.stress = self.cal_stress()
        
    def cal_strain_matix(self):
        for i in range(len(self.element)):
            nodes = self.nodes[self.element[i]]
            x1 = nodes[0,0]
            x2 = nodes[1,0]
            x3 = nodes[2,0]
            y1 = nodes[0,1]
            y2 = nodes[1,1]
            y3 = nodes[2,1]
            B = np.array([[y2-y3,0,y3-y1,0,y1-y2,0],
                          [0,x3-x2,0,x1-x3,0,x2-x1],
                          [x3-x2,y2-y3,x1-x3,y3-y1,x2-x1,y1-y2]])
            B = B/(2*self.cal_area(nodes))
            # print(B)
        return B
    
    def cal_area(self,nodes):
        x1 = nodes[0,0]
        x2 = nodes[1,0]
        x3 = nodes[2,0]
        y1 = nodes[0,1]
        y2 = nodes[1,1]
        y3 = nodes[2,1]
        area = 0.5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
        return area
    
    def cal_D(self):
        D = np.array([[1,self.nu,0],
                      [self.nu,1,0],
                      [0,0,(1-self.nu)/2]])
        D = self.E/(1-self.nu**2)*D
        return D
            
    def cal_stress(self):
        stress = []
        for i in range(len(self.element)):
            nodes = self.nodes[self.element[i]]
            # B = self.cal_strain_matix(nodes)
            stress.append(np.dot(self.D,np.dot(self.B,self.d)))
        return stress

