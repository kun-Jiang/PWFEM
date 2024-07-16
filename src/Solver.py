import numpy as np
from src.Boundary import Boundary
from src.Assemble import assemble


class solver:

    def __init__(self, steps, tolerance=1e-6, iter_max=10):
        
        self.K = assemble.K
        self.F = Boundary.Force
        self.steps = steps
        self.Penalty = Boundary.Penalty
        self.tolerance = tolerance
        self.iter_max = iter_max
        # for i in range(steps):
        self.solve()
        
            
    def solve(self):
        self.disp = np.zeros(len(self.F))
        for step in range(self.steps):
            # Initialize the load increment
            delta_F = self.F / self.steps
            residual = delta_F
            iter_count = 0
            # Iterative solution
            while np.linalg.norm(residual) > self.tolerance and iter_count < self.iter_max:
                # Solve for the incremental displacement
                delta_u = np.linalg.solve(self.K, residual)
                self.disp += delta_u
                # Update the residual
                internal_force = np.dot(self.K, delta_u)
                residual = (delta_F - internal_force)/Boundary.Penalty
                iter_count += 1
            if iter_count == self.iter_max:
                print(f"Warning: Load step {step+1} did not converge")
            
        return self.disp
        