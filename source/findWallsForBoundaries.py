
# This script is called by getBoundaryConditions.m, it takes flowRegion as input and outputs a list of walls and corners in the domain
# Each wall is defined by six columns in a list. The columns stand for [xi xf yi yf zi zf]
# Corners are defined by their indices (corners.limits) and they directions (corners.dir). For example, a corner defined by direction [-1 1 0] faces backwards in x, upwards in y and is constant in z
# corners.adiabatic indicates if the corner is adiabatic or isothermal. It defaults to isothermal.
# Regions that are completely confined inside walls are placed in the insideWalls variable

import numpy as np

# transformei este script em uma função (gigiaero - 29/08/2024)

def findWallsForBoundaries(flowRegion,mesh):

    # Find and remove infinitelly thin walls
    # in x
    flowRegion[:, 1:-1, :] |= (flowRegion[:, :-2, :] & flowRegion[:, 2:, :])

    # in y
    flowRegion[:, :, 1:-1] |= (flowRegion[:, :, :-2] & flowRegion[:, :, 2:])

    # in z
    flowRegion[1:-1, :, :] |= (flowRegion[:-2, :, :] & flowRegion[2:, :, :])

    # Find walls
    wallFront = np.diff(flowRegion, axis=1) == 1
    wallBack = np.diff(flowRegion, axis=1) == -1

    wallUp = np.diff(flowRegion, axis=2) == 1
    wallDown = np.diff(flowRegion, axis=2) == -1

    wallRight = np.diff(flowRegion, axis=3) == 1
    wallLeft = np.diff(flowRegion, axis=3) == -1

    # Find regions for boundaries in walls
    # In x
    wallFrontLimits = []
    wallBackLimits = []

    for i in range(mesh.nx):
        localWall = wallFront[i, :, :]  # Find front facing walls
        if np.any(localWall):
            for k in range(mesh.nz):
                if k == 0 or not wallStarts:  # Get starting and ending points in y
                    wallStarts = np.where([localWall[0, 0, k], np.diff(localWall[0, :, k], axis=1) == 1])[0]
                    wallEnds = np.where([np.diff(localWall[0, :, k], axis=1) == -1, localWall[0, -1, k]])[0]
                    kStart = k
                if k == mesh.nz - 1 or np.any(localWall[0, :, k] != localWall[0, :, k + 1]):  # If reached end in z or wall change, record boundary
                    for j in range(len(wallStarts)):
                        wallFrontLimits.append([i, i, wallStarts[j], wallEnds[j], kStart, kStart + 1])
                    wallStarts = []
                    wallEnds = []

        localWall = wallBack[i, :, :]  # Find down facing walls
        if np.any(localWall):
            for k in range(mesh.nz):
                if k == 0 or not wallStarts:  # Get starting and ending points in y
                    wallStarts = np.where([localWall[0, 0, k], np.diff(localWall[0, :, k], axis=1) == 1])[0]
                    wallEnds = np.where([np.diff(localWall[0, :, k], axis=1) == -1, localWall[0, -1, k]])[0]
                    kStart = k
                if k == mesh.nz - 1 or np.any(localWall[0, :, k] != localWall[0, :, k + 1]):  # If reached end in z or wall change, record boundary
                    for j in range(len(wallStarts)):
                        wallBackLimits.append([i, i, wallStarts[j], wallEnds[j], kStart, kStart + 1])
                    wallStarts = []
                    wallEnds = []

    # in y
    wallUpLimits = []
    wallDownLimits = []

    for j in range(mesh.ny):
        localWall = wallUp[:, j, :]  # Find front facing walls
        if np.any(localWall):
            for k in range(mesh.nz):
                if k == 0 or not wallStarts:  # Get starting and ending points in x
                    wallStarts = np.where([localWall[0, 0, k], np.diff(localWall[:, 0, k], axis=0) == 1])[0]
                    wallEnds = np.where([np.diff(localWall[:, 0, k], axis=0) == -1, localWall[-1, 0, k]])[0]
                    kStart = k
                if k == mesh.nz - 1 or np.any(localWall[:, 0, k] != localWall[:, 0, k + 1]):  # If reached end in z or wall change, record boundary
                    for i in range(len(wallStarts)):
                        wallUpLimits.append([wallStarts[i], wallEnds[i], j, j, kStart, kStart + 1])
                    wallStarts = []
                    wallEnds = []

        localWall = wallDown[:, j, :]  # Find front facing walls
        if np.any(localWall):
            for k in range(mesh.nz):
                if k == 0 or not wallStarts:  # Get starting and ending points in x
                    wallStarts = np.where([localWall[0, 0, k], np.diff(localWall[:, 0, k], axis=0) == 1])[0]
                    wallEnds = np.where([np.diff(localWall[:, 0, k], axis=0) == -1, localWall[-1, 0, k]])[0]
                    kStart = k
                if k == mesh.nz - 1 or np.any(localWall[:, 0, k] != localWall[:, 0, k + 1]):  # If reached end in z or wall change, record boundary
                    for i in range(len(wallStarts)):
                        wallDownLimits.append([wallStarts[i], wallEnds[i], j, j, kStart, kStart + 1])
                    wallStarts = []
                    wallEnds = []

    # in z
    wallRightLimits = []
    wallLeftLimits = []

    for k in range(mesh.nz):
        localWall = wallRight[:, :, k]  # Find left facing walls
        if np.any(localWall):
            for i in range(mesh.nx):
                if i == 0 or not wallStarts:  # Get starting and ending points in z
                    wallStarts = np.where([localWall[i, 0, 0], np.diff(localWall[i, :, 0], axis=1) == 1])[0]
                    wallEnds = np.where([np.diff(localWall[i, :, 0], axis=1) == -1, localWall[i, -1, 0]])[0]
                    iStart = i
                if i == mesh.nx - 1 or np.any(localWall[i, :, 0] != localWall[i + 1, :, 0]):  # If reached end in x or wall change, record boundary
                    for j in range(len(wallStarts)):
                        wallRightLimits.append([iStart, iStart + 1, wallStarts[j], wallEnds[j], k, k + 1])
                    wallStarts = []
                    wallEnds = []

        localWall = wallLeft[:, :, k]  # Find right facing walls
        if np.any(localWall):
            for i in range(mesh.nx):
                if i == 0 or not wallStarts:  # Get starting and ending points in z
                    wallStarts = np.where([localWall[i, 0, 0], np.diff(localWall[i, :, 0], axis=1) == 1])[0]
                    wallEnds = np.where([np.diff(localWall[i, :, 0], axis=1) == -1, localWall[i, -1, 0]])[0]
                    iStart = i
                if i == mesh.nx - 1 or np.any(localWall[i, :, 0] != localWall[i + 1, :, 0]):  # If reached end in x or wall change, record boundary
                    for j in range(len(wallStarts)):
                        wallLeftLimits.append([iStart, iStart + 1, wallStarts[j], wallEnds[j], k, k + 1])
                    wallStarts = []
                    wallEnds = []

    # Join adjacent walls
    wallFrontLimits = np.array(wallFrontLimits)
    wallBackLimits = np.array(wallBackLimits)
    wallUpLimits = np.array(wallUpLimits)
    wallDownLimits = np.array(wallDownLimits)

    for wallLimits in [wallFrontLimits, wallBackLimits, wallUpLimits, wallDownLimits]:
        done = False
        while not done:
            done = True
            nWalls = wallLimits.shape[0]
            for i in range(nWalls - 1):
                for j in range(i + 1, nWalls):
                    if np.all(wallLimits[i, :4] == wallLimits[j, :4]) and wallLimits[i, 5] == wallLimits[j, 4] - 1:
                        toMerge = np.array([i, j])
                        done = False
                        break
                if not done:
                    break
            if not done:
                wallLimits[toMerge[0], 5] = wallLimits[toMerge[1], 5]
                wallLimits = np.delete(wallLimits, toMerge[1], axis=0)

    # Find corners
    corners = {"limits": [], "dir": [], "adiabatic": np.zeros(len(wallFrontLimits), dtype=bool)}

    # With 2 walls
    # Constant z
    cornersMatrix = np.zeros((mesh.nx, mesh.ny, mesh.nz, 4), dtype=bool)
    cornersMatrix[:, :, :, 0] = wallFront & wallUp
    cornersMatrix[:, :, :, 1] = wallFront & wallDown
    cornersMatrix[:, :, :, 2] = wallBack & wallUp
    cornersMatrix[:, :, :, 3] = wallBack & wallDown
    cornerDirections = np.array([[1, 1, 0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0]])

    for i in range(mesh.nx):
        for j in range(mesh.ny):
            for m in range(4):
                cornerRow = cornersMatrix[i, j, :, m]
                if np.any(cornerRow):
                    cornerStarts = np.where(np.diff(cornerRow) == 1)[0] + 1
                    cornerEnds = np.where(np.diff(cornerRow) == -1)[0]
                    for n in range(len(cornerStarts)):
                        corners["limits"].append([i, i, j, j, cornerStarts[n], cornerEnds[n]])
                        corners["dir"].append(cornerDirections[m])

    # Constant x
    cornersMatrix = np.zeros((mesh.nx, mesh.ny, mesh.nz, 4), dtype=bool)
    cornersMatrix[:, :, :, 0] = wallFront & wallRight
    cornersMatrix[:, :, :, 1] = wallFront & wallLeft
    cornersMatrix[:, :, :, 2] = wallBack & wallRight
    cornersMatrix[:, :, :, 3] = wallBack & wallLeft
    cornerDirections = np.array([[1, 0, 1], [1, 0, -1], [-1, 0, 1], [-1, 0, -1]])

    for i in range(mesh.nx):
        for k in range(mesh.nz):
            for m in range(4):
                cornerRow = cornersMatrix[i, :, k, m]
                if np.any(cornerRow):
                    cornerStarts = np.where(np.diff(cornerRow) == 1)[0] + 1
                    cornerEnds = np.where(np.diff(cornerRow) == -1)[0]
                    for n in range(len(cornerStarts)):
                        corners["limits"].append([i, i, cornerStarts[n], cornerEnds[n], k, k])
                        corners["dir"].append(cornerDirections[m])

    # Constant y
    cornersMatrix = np.zeros((mesh.nx, mesh.ny, mesh.nz, 4), dtype=bool)
    cornersMatrix[:, :, :, 0] = wallUp & wallRight
    cornersMatrix[:, :, :, 1] = wallUp & wallLeft
    cornersMatrix[:, :, :, 2] = wallDown & wallRight
    cornersMatrix[:, :, :, 3] = wallDown & wallLeft
    cornerDirections = np.array([[0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1]])

    for j in range(mesh.ny):
        for k in range(mesh.nz):
            for m in range(4):
                cornerRow = cornersMatrix[:, j, k, m]
                if np.any(cornerRow):
                    cornerStarts = np.where(np.diff(cornerRow) == 1)[0] + 1
                    cornerEnds = np.where(np.diff(cornerRow) == -1)[0]
                    for n in range(len(cornerStarts)):
                        corners["limits"].append([cornerStarts[n], cornerEnds[n], j, j, k, k])
                        corners["dir"].append(cornerDirections[m])

    # With 3 walls
    cornersMatrix = np.zeros((mesh.nx, mesh.ny, mesh.nz, 8), dtype=bool)
    cornersMatrix[:, :, :, 0] = wallFront & wallUp & wallRight
    cornersMatrix[:, :, :, 1] = wallFront & wallUp & wallLeft
    cornersMatrix[:, :, :, 2] = wallFront & wallDown & wallRight
    cornersMatrix[:, :, :, 3] = wallFront & wallDown & wallLeft
    cornersMatrix[:, :, :, 4] = wallBack & wallUp & wallRight
    cornersMatrix[:, :, :, 5] = wallBack & wallUp & wallLeft
    cornersMatrix[:, :, :, 6] = wallBack & wallDown & wallRight
    cornersMatrix[:, :, :, 7] = wallBack & wallDown & wallLeft

    cornerDirections = np.array([[1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1], [-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]])

    I, J, K, M = np.where(cornersMatrix)
    corners["limits"].extend(list(zip(I, I, J, J, K, K)))
    corners["dir"].extend(cornerDirections[M])

    # Find regions completely contained inside walls
    isWall = ~flowRegion

    # Find wall limits in x
    insideWalls = []
    for k in range(mesh.nz):
        for j in range(mesh.ny):
            wallBegins = np.where([isWall[0, j, k], np.diff(isWall[:, j, k]) == 1])[0]
            wallEnds = np.where([np.diff(isWall[:, j, k]) == -1, isWall[-1, j, k]])[0]
            for i in range(len(wallBegins)):
                insideWalls.append([wallBegins[i], wallEnds[i], j, j, k, k])

    # Sort by same wall limits
    if len(insideWalls) > 0:
        _, ind = np.unique(np.array(insideWalls)[:, :4], axis=0, return_index=True, return_inverse=False)
        ind = np.argsort(ind)
        insideWalls = np.array(insideWalls)[ind]

    # Merge limits in y
    i = 0
    while i < len(insideWalls):
        if (
            np.all(insideWalls[i, :4] == insideWalls[i + 1, :4])
            and insideWalls[i, 3] + 1 == insideWalls[i + 1, 2]
        ):
            insideWalls[i, 3] = insideWalls[i + 1, 3]
            insideWalls = np.delete(insideWalls, i + 1, axis=0)
        else:
            i += 1

    # Merge limits in z
    i = 0
    while i < len(insideWalls):
        if (
            np.all(insideWalls[i, :4] == insideWalls[i + 1, :4])
            and insideWalls[i, 4] + 1 == insideWalls[i + 1, 3]
        ):
            insideWalls[i, 4] = insideWalls[i + 1, 4]
            insideWalls = np.delete(insideWalls, i + 1, axis=0)
        else:
            i += 1

    return corners,insideWalls,wallFrontLimits,wallBackLimits,wallUpLimits,wallDownLimits,wallRightLimits,wallLeftLimits