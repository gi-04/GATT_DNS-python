import numpy as np

# This script is called by makeMatrices.m and gets the flowRegion array as input and outputs three arrays, one for each direction
# For example, typeMapX is ny by nz and each entry contains the type of the derivative to be used in the row
# derivStartsX and derivEndsX contain the beginning and ending of each derivative type in typeMapX

## In X
# Get matrix
flowRegion = np.transpose(boundary.flowRegion, (1, 0, 2))

# Allocate variables
derivStartsX = []
derivEndsX = []

# Find all different types of geometries
C = []
for k in range(mesh.nz):
    C.append(np.unique(flowRegion[:, :, k], axis=0))
C = np.unique(np.vstack(C), axis=0)

# Identify where these regions are present
typeMapX = np.zeros((mesh.ny, mesh.nz), dtype=int)
for j in range(mesh.ny):
    for k in range(mesh.nz):
        typeMapX[j, k] = np.where(np.all(flowRegion[j, :, k] == C, axis=1))[0][0] + 1

for i in range(C.shape[0]):  # Find walls in each of the regions types
    if np.any(C[i, :] != 0):  # Only record regions with flow
        if mesh.x.periodic:
            derivStartsX.append(np.unique(np.where(np.diff(C[i, :]) == 1)[0]))
            derivEndsX.append(np.unique(np.where(np.diff(C[i, :]) == -1)[0] + 1))
        else:
            derivStartsX.append(np.unique(np.concatenate(([0], np.where(np.diff(C[i, :]) == 1)[0]))))
            derivEndsX.append(np.unique(np.concatenate((np.where(np.diff(C[i, :]) == -1)[0] + 1, [mesh.nx]))))
    else:
        typeMapX[typeMapX == i + 1] = 0
        typeMapX[typeMapX > i + 1] -= 1

typeMapX[typeMapX == 0] = np.max(typeMapX)

if hasattr(mesh.x, 'breakPoint'):
    for j in range(mesh.x.breakPoint.shape[0]):
        ind = mesh.x.breakPoint[j, 1:5]
        types = np.unique(typeMapX[ind[0]:ind[1], ind[2]:ind[3]])
        for t in types:
            typeN = np.max(typeMapX) + 1
            derivStartsX.append(derivStartsX[t - 1])
            derivEndsX.append(derivEndsX[t - 1])

            derivStartsX[typeN - 1] = np.append(derivStartsX[typeN - 1], mesh.x.breakPoint[j, 0] + 1)
            derivEndsX[typeN - 1] = np.append(derivEndsX[typeN - 1], mesh.x.breakPoint[j, 0])

            typeMapTemp = typeMapX[ind[0]:ind[1], ind[2]:ind[3]]
            typeMapTemp[typeMapTemp == t] = typeN
            typeMapX[ind[0]:ind[1], ind[2]:ind[3]] = typeMapTemp

## In Y
# Get matrix
flowRegion = boundary.flowRegion

# Allocate variables
derivStartsY = []
derivEndsY = []

# Find all different types of geometries
C = []
for k in range(mesh.nz):
    C.append(np.unique(flowRegion[:, :, k], axis=0))
C = np.unique(np.vstack(C), axis=0)

# Identify where these regions are present
typeMapY = np.zeros((mesh.nx, mesh.nz), dtype=int)
for i in range(mesh.nx):
    for k in range(mesh.nz):
        typeMapY[i, k] = np.where(np.all(flowRegion[i, :, k] == C, axis=1))[0][0] + 1

for i in range(C.shape[0]):  # Find walls in each of the regions types
    if np.any(C[i, :] != 0):  # Only record regions with flow
        if mesh.y.periodic:
            derivStartsY.append(np.unique(np.where(np.diff(C[i, :]) == 1)[0]))
            derivEndsY.append(np.unique(np.where(np.diff(C[i, :]) == -1)[0] + 1))
        else:
            derivStartsY.append(np.unique(np.concatenate(([0], np.where(np.diff(C[i, :]) == 1)[0]))))
            derivEndsY.append(np.unique(np.concatenate((np.where(np.diff(C[i, :]) == -1)[0] + 1, [mesh.ny]))))
    else:
        typeMapY[typeMapY == i + 1] = 0
        typeMapY[typeMapY > i + 1] -= 1

typeMapY[typeMapY == 0] = np.max(typeMapY)

if hasattr(mesh.y, 'breakPoint'):
    for j in range(mesh.y.breakPoint.shape[0]):
        ind = mesh.y.breakPoint[j, 1:5]
        types = np.unique(typeMapY[ind[0]:ind[1], ind[2]:ind[3]])
        for t in types:
            typeN = np.max(typeMapY) + 1
            derivStartsY.append(derivStartsY[t - 1])
            derivEndsY.append(derivEndsY[t - 1])

            derivStartsY[typeN - 1] = np.append(derivStartsY[typeN - 1], mesh.y.breakPoint[j, 0] + 1)
            derivEndsY[typeN - 1] = np.append(derivEndsY[typeN - 1], mesh.y.breakPoint[j, 0])

            typeMapTemp = typeMapY[ind[0]:ind[1], ind[2]:ind[3]]
            typeMapTemp[typeMapTemp == t] = typeN
            typeMapY[ind[0]:ind[1], ind[2]:ind[3]] = typeMapTemp

## In Z
# Get matrix
flowRegion = np.transpose(boundary.flowRegion, (0, 2, 1))

# Allocate variables
derivStartsZ = []
derivEndsZ = []

# Find all different types of geometries
C = []
for j in range(mesh.ny):
    C.append(np.unique(flowRegion[:, :, j], axis=0))
C = np.unique(np.vstack(C), axis=0)

# Identify where these regions are present
typeMapZ = np.zeros((mesh.nx, mesh.ny), dtype=int)
for i in range(mesh.nx):
    for j in range(mesh.ny):
        typeMapZ[i, j] = np.where(np.all(flowRegion[i, :, j] == C, axis=1))[0][0] + 1

for i in range(C.shape[0]):  # Find walls in each of the regions types
    if np.any(C[i, :] != 0):  # Only record regions with flow
        if mesh.z.periodic:
            derivStartsZ.append(np.unique(np.where(np.diff(C[i, :]) == 1)[0]))
            derivEndsZ.append(np.unique(np.where(np.diff(C[i, :]) == -1)[0] + 1))
        else:
            derivStartsZ.append(np.unique(np.concatenate(([0], np.where(np.diff(C[i, :]) == 1)[0]))))
            derivEndsZ.append(np.unique(np.concatenate((np.where(np.diff(C[i, :]) == -1)[0] + 1, [mesh.nz]))))
    else:
        typeMapZ[typeMapZ == i + 1] = 0
        typeMapZ[typeMapZ > i + 1] -= 1

typeMapZ[typeMapZ == 0] = np.max(typeMapZ)

if hasattr(mesh.z, 'breakPoint'):
    for j in range(mesh.z.breakPoint.shape[0]):
        ind = mesh.z.breakPoint[j, 1:5]
        types = np.unique(typeMapZ[ind[0]:ind[1], ind[2]:ind[3]])
        for t in types:
            typeN = np.max(typeMapZ) + 1
            derivStartsZ.append(derivStartsZ[t - 1])
            derivEndsZ.append(derivEndsZ[t - 1])

            derivStartsZ[typeN - 1] = np.append(derivStartsZ[typeN - 1], mesh.z.breakPoint[j, 0] + 1)
            derivEndsZ[typeN - 1] = np.append(derivEndsZ[typeN - 1], mesh.z.breakPoint[j, 0])

            typeMapTemp = typeMapZ[ind[0]:ind[1], ind[2]:ind[3]]
            typeMapTemp[typeMapTemp == t] = typeN
            typeMapZ[ind[0]:ind[1], ind[2]:ind[3]] = typeMapTemp



