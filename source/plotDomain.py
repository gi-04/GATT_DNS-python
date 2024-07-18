# This script plots the domain before runtime

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

colors = ['b', 'g', 'r', 'c', 'm', 'y']

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.hold(True)

for i in range(6):
    if i == 0:
        wallToPlot = boundary.wall.front
    elif i == 1:
        wallToPlot = boundary.wall.back
    elif i == 2:
        wallToPlot = boundary.wall.up
    elif i == 3:
        wallToPlot = boundary.wall.down
    elif i == 4:
        wallToPlot = boundary.wall.right
    else:
        wallToPlot = boundary.wall.left
    
    if wallToPlot is not None:
        if i in [0, 1]:
            X = mesh.X[wallToPlot[:, [0, 0, 0, 0]]]
            Y = mesh.Y[wallToPlot[:, [2, 3, 3, 2]]]
            Z = mesh.Z[wallToPlot[:, [4, 4, 5, 5]]]
        elif i in [2, 3]:
            X = mesh.X[wallToPlot[:, [0, 1, 1, 0]]]
            Y = mesh.Y[wallToPlot[:, [2, 2, 2, 2]]]
            Z = mesh.Z[wallToPlot[:, [4, 4, 5, 5]]]
        else:
            X = mesh.X[wallToPlot[:, [0, 1, 1, 0]]]
            Y = mesh.Y[wallToPlot[:, [2, 2, 3, 3]]]
            Z = mesh.Z[wallToPlot[:, [4, 4, 4, 4]]]
        
        for j in range(X.shape[0]):
            ax.patch(X[j, :], Z[j, :], Y[j, :], color=colors[i])

corners = boundary.corners.limits
cornerDir = boundary.corners.dir
for i in range(corners.shape[0]):
    if np.sum(np.abs(cornerDir[i, :])) == 2:
        ax.plot3D(mesh.X[corners[i, [0, 1]]], mesh.Z[corners[i, [4, 5]]], mesh.Y[corners[i, [2, 3]]], 'r', linewidth=2)
    else:
        ax.plot3D(mesh.X[corners[i, 0]], mesh.Z[corners[i, 4]], mesh.Y[corners[i, 2]], 'ro')
    
    X = np.mean(mesh.X[corners[i, [0, 1]]])
    Y = np.mean(mesh.Y[corners[i, [2, 3]]])
    Z = np.mean(mesh.Z[corners[i, [4, 5]]])
    
    scale = 0.1
    ax.plot3D([X, X + cornerDir[i, 0] * scale], [Z, Z + cornerDir[i, 2] * scale], [Y, Y + cornerDir[i, 1] * scale], 'g', linewidth=2)

ax.set_aspect('equal')
ax.view_init(30, 30)
plt.tight_layout()
plt.show()


