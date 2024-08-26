import numpy as np

# transformei este script em uma função. e notei que o nome deste arquivo é "calcSFDregion.py". será
# que deveria ser alterado pra calcSFDRegion.py pra seguir o padrão dos nomes dos demais arquivos?
#  (gigiaero - 26/08/2024)

def calcSFDRegion(mesh,numMethods):

    # This script is called if the SFD is to be applied only at the buffer zones.
    # It generates an array the same size as the domain which smoothly goes from 0 at the physical region to SFD_X at the buffer zones.

    SFD_X = np.ones((mesh.nx, mesh.ny, mesh.nz))

    # Set default value for applyX, applyY, applyZ
    if not hasattr(numMethods.SFD, 'applyX') or (isinstance(numMethods.SFD.applyX, bool) and numMethods.SFD.applyX):
        numMethods.SFD.applyX = 2
    if not hasattr(numMethods.SFD, 'applyY') or (isinstance(numMethods.SFD.applyY, bool) and numMethods.SFD.applyY):
        numMethods.SFD.applyY = 2
    if not hasattr(numMethods.SFD, 'applyZ') or (isinstance(numMethods.SFD.applyZ, bool) and numMethods.SFD.applyZ):
        numMethods.SFD.applyZ = 2

    if numMethods.SFD.applyX == 2 or numMethods.SFD.applyX == -1:
        for i in range(1, mesh.x.buffer.i.n + 1):
            SFD_X[i-1, :, :] *= (0.5 - 0.5 * np.cos(np.pi * (i-1) / mesh.x.buffer.i.n))

    if numMethods.SFD.applyX == 2 or numMethods.SFD.applyX == 1:
        for i in range(1, mesh.x.buffer.f.n + 1):
            SFD_X[-i, :, :] *= (0.5 - 0.5 * np.cos(np.pi * (i-1) / mesh.x.buffer.f.n))

    if numMethods.SFD.applyY == 2 or numMethods.SFD.applyY == -1:
        for j in range(1, mesh.y.buffer.i.n + 1):
            SFD_X[:, j-1, :] *= (0.5 - 0.5 * np.cos(np.pi * (j-1) / mesh.y.buffer.i.n))

    if numMethods.SFD.applyY == 2 or numMethods.SFD.applyY == 1:
        for j in range(1, mesh.y.buffer.f.n + 1):
            SFD_X[:, -j, :] *= (0.5 - 0.5 * np.cos(np.pi * (j-1) / mesh.y.buffer.f.n))

    if numMethods.SFD.applyZ == 2 or numMethods.SFD.applyZ == -1:
        for k in range(1, mesh.z.buffer.i.n + 1):
            SFD_X[:, :, k-1] *= (0.5 - 0.5 * np.cos(np.pi * (k-1) / mesh.z.buffer.i.n))

    if numMethods.SFD.applyZ == 2 or numMethods.SFD.applyZ == 1:
        for k in range(1, mesh.z.buffer.f.n + 1):
            SFD_X[:, :, -k] *= (0.5 - 0.5 * np.cos(np.pi * (k-1) / mesh.z.buffer.f.n))

    SFD_X = numMethods.SFD.X * (1 - SFD_X)

    if hasattr(numMethods.SFD, 'extraRegion'):
        for ER in numMethods.SFD.extraRegion:
            # Get radius from center point
            R = np.sqrt((mesh.X.T - ER.location[0])**2 / ER.size[0]**2 + 
                        (mesh.Y - ER.location[1])**2 / ER.size[1]**2 + 
                        (np.transpose(mesh.Z, (0, 2, 1)) - ER.location[2])**2 / ER.size[2]**2)
            R[R > 1] = 1
            R = 0.5 + 0.5 * np.cos(np.pi * R)
            SFD_X += ER.X * R
    
    return numMethods, SFD_X



