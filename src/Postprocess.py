import matplotlib.pyplot as plt
import matplotlib.tri as mtri

class PostProcess:
    def __init__(self,nodes,elements,result):
        self.nodes = nodes
        self.element = elements
        self.stress = result
        self.plot_stress()
    
    def plot_stress(self):
        x = self.nodes[:, 0]
        y = self.nodes[:, 1]
        z = self.stress[0]
        triang = mtri.Triangulation(self.nodes[:, 0], self.nodes[:, 1],self.element)
        fig = plt.figure()
        plt.triplot(triang,'k.-',lw=1)
        plt.tripcolor(triang,z,shading='gouraud')
        plt.tight_layout()
        plt.gca().set_aspect('equal')
        plt.set_cmap('jet')
        plt.colorbar(mappable=None, cax=None, ax=None, label="stress")
        plt.show()