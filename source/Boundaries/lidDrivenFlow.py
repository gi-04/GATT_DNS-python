import numpy as np

# Define flow region
# Find which nodes will actually contain a flow and which ones will be in or at a wall
flowRegion = np.ones((mesh.nx, mesh.ny, mesh.nz), dtype=bool)

# Add cavities to the flow region
if hasattr(flowType, 'cav'):
    for cavity in flowType.cav:
        x, y, z = cavity.x, cavity.y, cavity.z
        flowRegion[(mesh.X > x[0]) & (mesh.X < x[1]) & 
                   (mesh.Y > y[0]) & (mesh.Y < y[1]) & 
                   (mesh.Z > z[0]) & (mesh.Z < z[1])] = True

# Remove roughnesses from the flow
if hasattr(flowType, 'rug'):
    for roughness in flowType.rug:
        x, y, z = roughness.x, roughness.y, roughness.z
        flowRegion[(mesh.X >= x[0]) & (mesh.X <= x[1]) & 
                   (mesh.Y >= y[0]) & (mesh.Y <= y[1]) & 
                   (mesh.Z >= z[0]) & (mesh.Z <= z[1])] = False

# Add outer walls
flowRegion[0, :, :] = False
flowRegion[-1, :, :] = False
flowRegion[:, 0, :] = False
flowRegion[:, -1, :] = False

# Get walls
findWallsForBoundaries()

# Add walls to boundary conditions
wall_directions = ['xi', 'xf', 'yi', 'yf', 'zi', 'zf']
wall_positions = [wallFrontLimits, wallBackLimits, wallUpLimits, 
                  wallDownLimits, wallRightLimits, wallLeftLimits]

var, type_, dir_, val, xi, xf, yi, yf, zi, zf = [], [], [], [], [], [], [], [], [], []

for wallDir, wallPosition in zip(wall_directions, wall_positions):
    for wall in wallPosition:
        var.append('p')
        type_.append('neu')
        dir_.append(wallDir)
        val.append(0)
        xi.append(wall[0])
        xf.append(wall[1])
        yi.append(wall[2])
        yf.append(wall[3])
        zi.append(wall[4])
        zf.append(wall[5])

# Add regions that are inside walls
for wall in insideWalls:
    for variable in ['p', 'u', 'v', 'w']:
        var.append(variable)
        type_.append('dir')
        dir_.append('yi')
        val.append(P0 if variable == 'p' else 0)
        xi.append(wall[0])
        xf.append(wall[1])
        yi.append(wall[2])
        yf.append(wall[3])
        zi.append(wall[4])
        zf.append(wall[5])
    
    if hasattr(flowType, 'tWallRelative'):
        eWall = E0 * flowType.tWallRelative
    else:
        eWall = E0
    
    var.append('e')
    type_.append('dir')
    dir_.append('yi')
    val.append(eWall)
    xi.append(wall[0])
    xf.append(wall[1])
    yi.append(wall[2])
    yf.append(wall[3])
    zi.append(wall[4])
    zf.append(wall[5])


