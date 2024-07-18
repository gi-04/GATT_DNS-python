
import numpy as np

# Subroutine for an isothermal boundary layer with or without cavities and roughtnesses
# To be called from getBoundaryConditions

# Inflow
if mesh.X[0] <= 0:
    # u
    var.append('u')
    type.append('dir')
    dir.append('xi')
    val.append(1)
    xi.append(1)
    xf.append(1)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(1)
    zf.append(mesh.nz)

    # v
    var.append('v')
    type.append('dir')
    dir.append('xi')
    val.append(0)
    xi.append(1)
    xf.append(1)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(1)
    zf.append(mesh.nz)

    # w
    var.append('w')
    type.append('dir')
    dir.append('xi')
    val.append(0)
    xi.append(1)
    xf.append(1)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(1)
    zf.append(mesh.nz)

    # e
    var.append('e')
    type.append('dir')
    dir.append('xi')
    val.append(E0)
    xi.append(1)
    xf.append(1)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(1)
    zf.append(mesh.nz)
else:
    if not hasattr(flowType, 'disturb'):
        flowType.disturb = []
    else:
        flowType.disturb = flowType.disturb[1:] + [[]]
    flowType.disturb[0].x = [mesh.X[0], mesh.X[0]]
    flowType.disturb[0].y = [-np.inf, np.inf]
    flowType.disturb[0].z = [-np.inf, np.inf]
    flowType.disturb[0].var = 'UVRWE'
    flowType.disturb[0].type = 'holdInlet'
    flowType.disturb[0].active = True
    flowType.disturb[0].par = [0]  # Hold density

# p
var.append('p')
type.append('neu')
dir.append('xi')
val.append(0)
xi.append(1)
xf.append(1)
yi.append(1)
yf.append(mesh.ny)
zi.append(1)
zf.append(mesh.nz)

# Outflow
# u
var.append('u')
type.append('sec')
dir.append('xf')
val.append(0)
xi.append(mesh.nx)
xf.append(mesh.nx)
yi.append(1)
yf.append(mesh.ny)
zi.append(1)
zf.append(mesh.nz)

# v
var.append('v')
type.append('sec')
dir.append('xf')
val.append(0)
xi.append(mesh.nx)
xf.append(mesh.nx)
yi.append(1)
yf.append(mesh.ny)
zi.append(1)
zf.append(mesh.nz)

# w
var.append('w')
type.append('sec')
dir.append('xf')
val.append(0)
xi.append(mesh.nx)
xf.append(mesh.nx)
yi.append(1)
yf.append(mesh.ny)
zi.append(1)
zf.append(mesh.nz)

# e
var.append('e')
type.append('sec')
dir.append('xf')
val.append(0)
xi.append(mesh.nx)
xf.append(mesh.nx)
yi.append(1)
yf.append(mesh.ny)
zi.append(1)
zf.append(mesh.nz)

# p
var.append('p')
type.append('dir')
dir.append('xf')
val.append(P0)
xi.append(mesh.nx)
xf.append(mesh.nx)
yi.append(1)
yf.append(mesh.ny)
zi.append(1)
zf.append(mesh.nz)

# Outerflow
# u
var.append('u')
type.append('sec')
dir.append('yf')
val.append(0)
xi.append(1)
xf.append(mesh.nx)
yi.append(mesh.ny)
yf.append(mesh.ny)
zi.append(1)
zf.append(mesh.nz)

# v
var.append('v')
type.append('sec')
dir.append('yf')
val.append(0)
xi.append(1)
xf.append(mesh.nx)
yi.append(mesh.ny)
yf.append(mesh.ny)
zi.append(1)
zf.append(mesh.nz)

# w
var.append('w')
type.append('sec')
dir.append('yf')
val.append(0)
xi.append(1)
xf.append(mesh.nx)
yi.append(mesh.ny)
yf.append(mesh.ny)
zi.append(1)
zf.append(mesh.nz)

# e
var.append('e')
type.append('sec')
dir.append('yf')
val.append(0)
xi.append(1)
xf.append(mesh.nx)
yi.append(mesh.ny)
yf.append(mesh.ny)
zi.append(1)
zf.append(mesh.nz)

# p
var.append('p')
type.append('sec')
dir.append('yf')
val.append(0)
xi.append(1)
xf.append(mesh.nx)
yi.append(mesh.ny)
yf.append(mesh.ny)
zi.append(1)
zf.append(mesh.nz)

# Outerflow for non-periodic 3D
if mesh.nz > 1 and not mesh.z.periodic:
    # u
    var.append('u')
    type.append('neu')
    dir.append('zi')
    val.append(0)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(1)
    zf.append(1)

    # v
    var.append('v')
    type.append('neu')
    dir.append('zi')
    val.append(0)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(1)
    zf.append(1)

    # w
    var.append('w')
    type.append('neu')
    dir.append('zi')
    val.append(0)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(1)
    zf.append(1)

    # e
    var.append('e')
    type.append('neu')
    dir.append('zi')
    val.append(0)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(1)
    zf.append(1)

    # p
    var.append('p')
    type.append('neu')
    dir.append('zi')
    val.append(0)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(1)
    zf.append(1)
    
    # u
    var.append('u')
    type.append('neu')
    dir.append('zf')
    val.append(0)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(mesh.nz)
    zf.append(mesh.nz)

    # v
    var.append('v')
    type.append('neu')
    dir.append('zf')
    val.append(0)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(mesh.nz)
    zf.append(mesh.nz)

    # w
    var.append('w')
    type.append('neu')
    dir.append('zf')
    val.append(0)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(mesh.nz)
    zf.append(mesh.nz)

    # e
    var.append('e')
    type.append('neu')
    dir.append('zf')
    val.append(0)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(mesh.nz)
    zf.append(mesh.nz)

    # p
    var.append('p')
    type.append('neu')
    dir.append('zf')
    val.append(0)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(1)
    yf.append(mesh.ny)
    zi.append(mesh.nz)
    zf.append(mesh.nz)

# Define flow region
# Find which nodes will actually contain a flow and which ones will be in or at a wall
flowRegion = np.ones((mesh.nx, mesh.ny, mesh.nz), dtype=bool)

# Add flat plate
[~, wallJ] = np.min(np.abs(mesh.Y), axis=0)
flowRegion[:, :wallJ, :] = False

# Add cavities to the flow region
if hasattr(flowType, 'cav'):
    for i in range(len(flowType.cav)):
        x = flowType.cav[i].x
        y = flowType.cav[i].y
        z = flowType.cav[i].z
        
        flowRegion[mesh.X > x[0] & mesh.X < x[1], 
                  mesh.Y > y[0] & mesh.Y < y[1], 
                  mesh.Z > z[0] & mesh.Z < z[1]] = True

# Remove roughnesses from the flow
if hasattr(flowType, 'rug'):
    for i in range(len(flowType.rug)):
        x = flowType.rug[i].x
        y = flowType.rug[i].y
        z = flowType.rug[i].z
        
        flowRegion[mesh.X >= x[0] & mesh.X <= x[1], 
                  mesh.Y >= y[0] & mesh.Y <= y[1], 
                  mesh.Z >= z[0] & mesh.Z <= z[1]] = False

# Get walls
findWallsForBoundaries()

corners.adiabatic = np.ones(corners.dir.shape[0], dtype=int)

# Create wall for free-slip region
if mesh.X[0] < 0:
    wallUpLimits[:, 0] = np.maximum(wallUpLimits[:, 0], 
                                    np.where(mesh.X >= 0)[0][0])  # Move existing walls
    wallUpLimits = np.vstack((wallUpLimits, 
                              np.array([[1, mesh.X[mesh.X < 0][-1], 
                                        1, 1, 1, mesh.nz]])))
    
    mesh.x.breakPoint = np.array([mesh.X[mesh.X < 0][-1], 
                                 1, 1, 1, mesh.nz])


