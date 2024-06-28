import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
sys.path.append('/home/senseye3/study/c_math/pde/thermal_conductivity2D/thermal_conductivity2d')

from meshes.mesh import Point2D, Grid
from solver.solver import thermal_conductivity2D


corner_points = np.array([Point2D(0.,0.), Point2D(0.,0.001),
                         Point2D(0.01,0.001), Point2D(0.01,0.)])
grid = Grid(corner_points, 100, 10)
grid.fill_grid()

t = thermal_conductivity2D(grid, 100.)
t.set_init_condit(300.)
#t.set_boundary_condit(300., 300.)
t.solve()

# Generate x and y coordinates
x = np.linspace(0, 0.01, 100)
y = np.linspace(0, 0.001, 10)
X, Y = np.meshgrid(x, y)

# Generate random temperature data with matching dimensions to X and Y
u2D = np.zeros(( t.grid.M,t.grid.K))
for k in range(t.grid.K):
    for m in range(t.grid.M):
        i = m + k*t.grid.M
        u2D[m][k] = np.copy(t.u[i])
print(u2D)

# Create a heatmap of the temperature field with explicit X and Y coordinates
cs = plt.contourf(x, y, u2D, levels=30,
                  cmap=matplotlib.cm.magma_r)
cbar = plt.colorbar(cs)


# Add labels and title
plt.xlabel('X')
plt.ylabel('Y')
plt.title('2D Temperature Field with Explicit Coordinates')

# Show the plot
plt.show()