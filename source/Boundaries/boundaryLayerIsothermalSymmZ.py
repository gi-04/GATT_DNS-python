import numpy as np
from findWallsForBoundaries import findWallsForBoundaries

def boundaryLayerIsothermalSymmZ(mesh,var,type,dir,val,xi,xf,yi,yf,zi,zf,flowType,E0,P0):

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
        if 'disturb' not in flowType:
            flowType['disturb'] = []
        else:
            flowType['disturb'] = flowType['disturb'][1:]  # Shift elements

        flowType['disturb'].insert(0, {})
        flowType['disturb'][0]['x'] = [mesh.X[0], mesh.X[0]]
        flowType['disturb'][0]['y'] = [-np.inf, np.inf]
        flowType['disturb'][0]['z'] = [-np.inf, np.inf]
        flowType['disturb'][0]['var'] = 'UVRWE'
        flowType['disturb'][0]['type'] = 'holdInlet'
        flowType['disturb'][0]['active'] = True
        flowType['disturb'][0]['par'] = [0]

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

    # Define flow region
    # Find which nodes will actually contain a flow and which ones will be in or at a wall
    flowRegion = np.ones((mesh.nx, mesh.ny, mesh.nz), dtype=bool)

    # Add flat plate
    _, wallJ = np.min(np.abs(mesh.Y)), np.argmin(np.abs(mesh.Y))
    flowRegion[:, :wallJ, :] = False

    # Add cavities to the flow region
    if 'cav' in flowType:
        for i in range(len(flowType['cav'])):
            x = flowType['cav'][i]['x']
            y = flowType['cav'][i]['y']
            z = flowType['cav'][i]['z']
            
            flowRegion[(mesh.X > x[0]) & (mesh.X < x[1]), 
                        (mesh.Y > y[0]) & (mesh.Y < y[1]), 
                        (mesh.Z > z[0]) & (mesh.Z < z[1])] = True

    # Remove roughnesses from the flow
    if 'rug' in flowType:
        for i in range(len(flowType['rug'])):
            x = flowType['rug'][i]['x']
            y = flowType['rug'][i]['y']
            z = flowType['rug'][i]['z']
            
            flowRegion[(mesh.X >= x[0]) & (mesh.X <= x[1]), 
                        (mesh.Y >= y[0]) & (mesh.Y <= y[1]), 
                        (mesh.Z >= z[0]) & (mesh.Z <= z[1])] = False

    # Get walls
    corners,insideWalls,wallFrontLimits,wallBackLimits,wallUpLimits,wallDownLimits,wallRightLimits,wallLeftLimits = findWallsForBoundaries(flowRegion,mesh)

    # Create wall for free-slip region
    if mesh.X[0] < 0:
        wallUpLimits[:, 0] = np.maximum(wallUpLimits[:, 0], np.where(mesh.X >= 0)[0][0])  # Move existing walls
        wallUpLimits = np.vstack((wallUpLimits, [1, np.where(mesh.X < 0)[0][-1], 
                                                np.array([1, 1]) * np.where(mesh.Y >= 0)[0][0], 
                                                1, mesh.nz]))

        mesh.x.breakPoint = [np.where(mesh.X < 0)[0][-1], 
                            np.array([1, 1]) * np.where(mesh.Y >= 0)[0][0], 
                            1, mesh.nz]

    # Add walls to boundary conditions
    for i in range(1, 7):
        if i == 1:
            wallPosition = wallFrontLimits
            wallDir = 'xi'
        elif i == 2:
            wallPosition = wallBackLimits
            wallDir = 'xf'
        elif i == 3:
            wallPosition = wallUpLimits
            wallDir = 'yi'
        elif i == 4:
            wallPosition = wallDownLimits
            wallDir = 'yf'
        elif i == 5:
            wallPosition = wallRightLimits
            wallDir = 'zi'
        elif i == 6:
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

    return var,type,dir,val,xi,xf,yi,yf,zi,zf,flowType,flowRegion,corners,wallFrontLimits,wallBackLimits,wallUpLimits,wallDownLimits,wallRightLimits,wallLeftLimits
