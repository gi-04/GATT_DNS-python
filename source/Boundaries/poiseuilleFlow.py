import numpy as np
from findWallsForBoundaries import findWallsForBoundaries

def poiseulleFlow(mesh,var,type,dir,val,xi,xf,yi,yf,zi,zf,flowType,E0,P0,flowParameters):

    # Define flow region
    # Find which nodes will actually contain a flow and which ones will be in or at a wall
    flowRegion = np.ones((mesh.nx, mesh.ny, mesh.nz), dtype=bool)

    # Add cavities to the flow region
    if hasattr(flowType, 'cav'):
        for cav in flowType.cav:
            x = cav.x
            y = cav.y
            z = cav.z
            
            flowRegion[(mesh.X > x[0]) & (mesh.X < x[1]), 
                    (mesh.Y > y[0]) & (mesh.Y < y[1]), 
                    (mesh.Z > z[0]) & (mesh.Z < z[1])] = True

    # Remove roughnesses from the flow
    if hasattr(flowType, 'rug'):
        for rug in flowType.rug:
            x = rug.x
            y = rug.y
            z = rug.z
            
            flowRegion[(mesh.X >= x[0]) & (mesh.X <= x[2]), 
                    (mesh.Y >= y[0]) & (mesh.Y <= y[2]), 
                    (mesh.Z >= z[0]) & (mesh.Z <= z[2])] = False

    # Add outer walls
    flowRegion[:, [0, -1], :] = False

    # Get walls
    corners,insideWalls,wallFrontLimits,wallBackLimits,wallUpLimits,wallDownLimits,wallRightLimits,wallLeftLimits,insideWalls = findWallsForBoundaries(flowRegion,mesh)

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

        if hasattr(flowType, 'tWallRelative'):
            eWall = E0 * flowType.tWallRelative
        else:
            eWall = E0

        var.append('e')
        type.append('dir')
        dir.append('yi')
        val.append(eWall)
        xi.append(insideWalls[i, 0])
        xf.append(insideWalls[i, 1])
        yi.append(insideWalls[i, 2])
        yf.append(insideWalls[i, 3])
        zi.append(insideWalls[i, 4])
        zf.append(insideWalls[i, 5])

    # Add moving walls
    var.append('u')
    type.append('dir')
    dir.append('yi')
    val.append(flowParameters.lowerWallVelocity)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(1)
    yf.append(1)
    zi.append(1)
    zf.append(mesh.nz)

    var.append('u')
    type.append('dir')
    dir.append('yf')
    val.append(flowParameters.upperWallVelocity)
    xi.append(1)
    xf.append(mesh.nx)
    yi.append(mesh.ny)
    yf.append(mesh.ny)
    zi.append(1)
    zf.append(mesh.nz)

    return corners,wallFrontLimits,wallBackLimits,wallUpLimits,wallDownLimits,wallRightLimits,wallLeftLimits
