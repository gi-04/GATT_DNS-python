
import numpy as np
from findWallsForBoundaries import findWallsForBoundaries

def periodicBox(mesh):

    # Subroutine for a periodic box
    # To be called from getBoundaryConditions

    flowRegion = np.ones((mesh.nx, mesh.ny, mesh.nz), dtype=bool)

    corners,_,wallFrontLimits,wallBackLimits,wallUpLimits,wallDownLimits,wallRightLimits,wallLeftLimits = findWallsForBoundaries(flowRegion,mesh)

    return corners,wallFrontLimits,wallBackLimits,wallUpLimits,wallDownLimits,wallRightLimits,wallLeftLimits

