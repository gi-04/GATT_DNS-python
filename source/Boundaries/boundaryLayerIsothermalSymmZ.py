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
        flowType.disturb.extend([])
    flowType.disturb.append({})
    flowType.disturb[0]['x'] = [mesh.X[0], mesh.X[0]]
    flowType.disturb[0]['y'] = [-float('inf'), float('inf')]
    flowType.disturb[0]['z'] = [-float('inf'), float('inf')]
    flowType.disturb[0]['var'] = 'UVRWE'
    flowType.disturb[0]['type'] = 'holdInlet'
    flowType.disturb[0]['active'] = True
    flowType.disturb[0]['par'] = [0]

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

# Outerflow for symmetry
if mesh.nz > 1:
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
    type.append('dir')
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
    type.append('dir')
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


