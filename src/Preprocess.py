import matplotlib.pyplot as plt
import matplotlib.tri as mtri

class Preprocess:
    def __init__(self, mesh):
        self.nodes = mesh.points
        self.element = mesh.cells_dict["triangle"]
    
    def mesh_visual(self):
        x = self.nodes[:, 0]
        y = self.nodes[:, 1]
        triang = mtri.Triangulation(x, y, self.element)
        fig = plt.figure()
        # ax = fig.add_subplot(111, projection="3d")
        for index,_ in enumerate(self.nodes):
            plt.annotate(index,(x[index],y[index]),color="red")
        plt.triplot(triang,'k.-',lw=1)
        plt.gca().set_aspect('equal')
        plt.show()