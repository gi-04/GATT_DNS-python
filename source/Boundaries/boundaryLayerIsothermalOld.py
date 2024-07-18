# Subroutine for an isothermal boundary layer with or without cavities and roughnesses
# To be called from getBoundaryConditions

# Inflow
# u
var = ['u']
type = ['dir']
dir = ['xi']
val = [1]
xi = [1]
xf = [1]
yi = [1]
yf = [mesh.ny]
zi = [1]
zf = [mesh.nz]

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
_, wallJ = min((abs(mesh.Y), range(len(mesh.Y))))
flowRegion[:, :wallJ, :] = False

# Add cavities to the flow region
if hasattr(flowType, 'cav'):
    for cav in flowType.cav:
        x = cav['x']
        y = cav['y']
        z = cav['z']
        flowRegion[(mesh.X > x[0]) & (mesh.X < x[1]), (mesh.Y > y[0]) & (mesh.Y < y[1]), (mesh.Z > z[0]) & (mesh.Z < z[1])] = True

# Remove roughnesses from the flow
if hasattr(flowType, 'rug'):
    for rug in flowType.rug:
        x = rug['x']
        y = rug['y']
        z = rug['z']
        flowRegion[(mesh.X >= x[0]) & (mesh.X <= x[2]), (mesh.Y >= y[0]) & (mesh.Y <= y[2]), (mesh.Z >= z[0]) & (mesh.Z <= z[2])] = False

# Get walls
findWallsForBoundaries()

# Create wall for free-slip region
if mesh.X[0] < 0:
    wallUpLimits[:, 0] = np.maximum(wallUpLimits[:, 0], np.searchsorted(mesh.X, 0, side='left'))
    wallUpLimits = np.vstack([wallUpLimits, [1, np.searchsorted(mesh.X, 0, side='right') - 1, np.searchsorted(mesh.Y, 0, side='left'), 1, mesh.nz]])

# Add walls to boundary conditions
for i in range(6):
    if i == 0:
        wallPosition = wallFrontLimits
        wallDir = 'xi'
    elif i == 1:
        wallPosition = wallBackLimits
        wallDir = 'xf'
    elif i == 2:
        wallPosition = wallUpLimits
        wallDir = 'yi'
    elif i == 3:
        wallPosition = wallDownLimits
        wallDir = 'yf'
    elif i == 4:
        wallPosition = wallRightLimits
        wallDir = 'zi'
    elif i == 5:
        wallPosition = wallLeftLimits
        wallDir = 'zf'

    for j in range(wallPosition.shape[0]):
        var.append('p')
        type.append('neu')
        dir.append(wallDir)
        val.append(0)
        xi.append(wallPosition[j, 0])
        xf.append(wallPosition[j, 1])
        yi.append(wallPosition[j, 2])
        yf.append(wallPosition[j, 3])
        zi.append(wallPosition[j, 4])
        zf.append(wallPosition[j, 5])

# Add free slip wall, it is the last up facing wall
if mesh.X[0] < 0:
    var.append('u')
    type.append('neu')
    dir.append('yi')
    val.append(0)
    xi.append(wallUpLimits[-1, 0])
    xf.append(wallUpLimits[-1, 1])
    yi.append(wallUpLimits[-1, 2])
    yf.append(wallUpLimits[-1, 3])
    zi.append(wallUpLimits[-1, 4])
    zf.append(wallUpLimits[-1, 5])

    var.append('w')
    type.append('neu')
    dir.append('yi')
    val.append(0)
    xi.append(wallUpLimits[-1, 0])
    xf.append(wallUpLimits[-1, 1])
    yi.append(wallUpLimits[-1, 2])
    yf.append(wallUpLimits[-1, 3])
    zi.append(wallUpLimits[-1, 4])
    zf.append(wallUpLimits[-1, 5])

# Add regions that are inside walls
for i in range(insideWalls.shape[0]):
    var.append('p')
    type.append('dir')
    dir.append('yi')
    val.append(P0)
    xi.append(insideWalls[i, 0])
    xf.append(insideWalls[i, 1])
    yi.append(insideWalls[i, 2])
    yf.append(insideWalls[i, 3])
    zi.append(insideWalls[i, 4])
    zf.append(insideWalls[i, 5])

    var.append('u')
    type.append('dir')
    dir.append('yi')
    val.append(0)
    xi.append(insideWalls[i, 0])
    xf.append(insideWalls[i, 1])
    yi.append(insideWalls[i, 2])
    yf.append(insideWalls[i, 3])
    zi.append(insideWalls[i, 4])
    zf.append(insideWalls[i, 5])

    var.append('v')
    type.append('dir')
    dir.append('yi')
    val.append(0)
    xi.append(insideWalls[i, 0])
    xf.append(insideWalls[i, 1])
    yi.append(insideWalls[i, 2])
    yf.append(insideWalls[i, 3])
    zi.append(insideWalls[i, 4])
    zf.append(insideWalls[i, 5])

    var.append('w')
    type.append('dir')
    dir.append('yi')
    val.append(0)
    xi.append(insideWalls[i, 0])
    xf.append(insideWalls[i, 1])
    yi.append(insideWalls[i, 2])
    yf.append(insideWalls[i, 3])
    zi.append(insideWalls[i, 4])
    zf.append(insideWalls[i, 5])

    var.append('e')
    type.append('dir')
    dir.append('yi')
    val.append(E0)
    xi.append(insideWalls[i, 0])
    xf.append(insideWalls[i, 1])
    yi.append(insideWalls[i, 2])
    yf.append(insideWalls[i, 3])
    zi.append(insideWalls[i, 4])
    zf.append(insideWalls[i, 5])



