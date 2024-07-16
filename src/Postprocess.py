import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

class Plot2D:
    def __init__(self, mesh, data, deformed=False, wireframe=False):
        self.deformed = deformed
        self.wireframe = wireframe
        nodes = mesh.points
        elements = mesh.cells_dict['quad']
        if self.deformed:
            self.PlotDeformed(nodes, elements, data)
        else:
            self.PlotUnDeformed(nodes, elements, data)
        
    def PlotUnDeformed(self,nodes, elements, data):
        fig, ax = plt.subplots(constrained_layout=True)
        x, y = nodes[:, 0], nodes[:, 1]
        X, Y = np.meshgrid(x, y)
        Z = griddata((x, y), data, (X, Y), method='linear')
        plt.contourf(X,Y,Z,levels=20,vmin=np.min(data),vmax=np.max(data),cmap='jet')
        cb = plt.colorbar()
        cb.set_ticks(np.linspace(min(data), max(data), 10))
        self.plot_wireframe(nodes, elements,data)
        plt.show()
        
    def PlotDeformed(self,nodes, elements, data):
        fig, ax = plt.subplots(constrained_layout=True)
        if self.deformed:
            nodes_deformed = nodes + data.reshape(-1, 2)
        x, y = nodes_deformed[:, 0], nodes_deformed[:, 1]
        X, Y = np.meshgrid(x, y)
        Z = griddata((x, y), data, (X, Y), method='linear')
        plt.contourf(X,Y,Z,levels=20,vmin=np.min(data),vmax=np.max(data),cmap='jet')
        cb = plt.colorbar()
        cb.set_ticks(np.linspace(min(data), max(data), 10))
        self.plot_wireframe(nodes, elements,data)
        plt.show()
        
    def plot_wireframe(self,nodes, elements, data):
        if self.deformed:
            nodes_deformed = nodes + data.reshape(-1, 2)
            nodes = nodes_deformed
        if self.wireframe:
            for element in elements:
                quad = nodes[element]
                plt.plot(*zip(*np.append(quad, quad[0:1], axis=0)), color='k', linewidth=0.5)
