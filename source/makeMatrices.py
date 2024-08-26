# Required python packages
import numpy as np
import scipy.sparse as sp
# Additional functions and files
from findDerivativeRegions import findDerivativeRegions
from finiteDifferenceCoefficients import finiteDifferenceCoefficients
from spatialFilterCoefficients import spatialFilterCoefficients

# foi preciso inicializar o dicionário matrices. quanto à função findDerivativeRegions,
# consertei explicitando as entradas e saídas. também corrigi nomes de variáveis pra
# ficarem iguais às do código de matlab (gigiaero - 26/08/2024)

def makeMatrices(mesh, domain, boundary, numMethods):
    """
    This function creates the matrices that will be used to compute derivatives and filters
    One pair of LHS and RHS will be computed for each type of region in the domain
    """

    # Find uniform regions
    typeMapX,typeMapY,typeMapZ,derivStartsX,derivStartsY,derivStartsZ,derivEndsX,derivEndsY,derivEndsZ = findDerivativeRegions(mesh,boundary)

    # Get finite differences coefficients
    centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS = finiteDifferenceCoefficients(numMethods.spatialDerivs)
    centeredStencilLHSb, centeredStencilRHSb, decenteredStencilLHSb, decenteredStencilRHSb = finiteDifferenceCoefficients(numMethods.spatialDerivs_buffer)
    filterStencilLhs, filterStencilRhs, filterDecenteredStencilLHS, filterDecenteredStencilRHS = spatialFilterCoefficients(numMethods.spatial_filter_strength, numMethods.filter_borders)

    # If the Z direction is not unitary but too small to fit the filter stencil, change the stencil just for it
    if mesh.nz > 1 and mesh.nz < 2 * len(filterStencilRhs) - 1:
        filterStencilLHSz = np.ones(mesh.nz)
        filterStencilRHSz = np.ones(mesh.nz)
        filterDecenteredStencilLHSz = np.ones(mesh.nz)
        filterDecenteredStencilRHSz = np.ones(mesh.nz)
    else:
        filterStencilLHSz = filterStencilLhs
        filterStencilRHSz = filterStencilRhs
        filterDecenteredStencilLHSz = filterDecenteredStencilLHS
        filterDecenteredStencilRHSz = filterDecenteredStencilRHS

    # If the mesh in Z has exactly 4 nodes and has no boundaries, switch it to spectral mode
    if mesh.nz == 4 and len(derivStartsZ) == 1 and len(derivStartsZ[0]) == 0 and len(derivEndsZ[0]) == 0:
        centeredStencilLHSz = np.ones(mesh.nz)
        centeredStencilRHSz = np.array([0, np.pi / 4])
        centeredStencilLHSbz = centeredStencilLHSz
        centeredStencilRHSbz = centeredStencilRHSz
    else:
        centeredStencilLHSz = centeredStencilLHS
        centeredStencilRHSz = centeredStencilRHS
        centeredStencilLHSbz = centeredStencilLHSb
        centeredStencilRHSbz = centeredStencilRHSb

    # Make matrices
    LHSx, RHSx = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, derivStartsX, derivEndsX, mesh.nx, [])
    LHSy, RHSy = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, derivStartsY, derivEndsY, mesh.ny, [])
    LHSz, RHSz = makeMatricesEachDirection(centeredStencilLHSz, centeredStencilRHSz, decenteredStencilLHS, decenteredStencilRHS, derivStartsZ, derivEndsZ, mesh.nz, [])

    LHSxb, RHSxb = makeMatricesEachDirection(centeredStencilLHSb, centeredStencilRHSb, decenteredStencilLHSb, decenteredStencilRHSb, derivStartsX, derivEndsX, mesh.nx, mesh.x.buffer)
    LHSyb, RHSyb = makeMatricesEachDirection(centeredStencilLHSb, centeredStencilRHSb, decenteredStencilLHSb, decenteredStencilRHSb, derivStartsY, derivEndsY, mesh.ny, mesh.y.buffer)
    LHSzb, RHSzb = makeMatricesEachDirection(centeredStencilLHSbz, centeredStencilRHSbz, decenteredStencilLHSb, decenteredStencilRHSb, derivStartsZ, derivEndsZ, mesh.nz, mesh.z.buffer)

    f_LHSx, f_RHSx = makeMatricesEachDirection(filterStencilLhs, filterStencilRhs, filterDecenteredStencilLHS, filterDecenteredStencilRHS, derivStartsX, derivEndsX, mesh.nx, [])
    f_LHSy, f_RHSy = makeMatricesEachDirection(filterStencilLhs, filterStencilRhs, filterDecenteredStencilLHS, filterDecenteredStencilRHS, derivStartsY, derivEndsY, mesh.ny, [])
    f_LHSz, f_RHSz = makeMatricesEachDirection(filterStencilLHSz, filterStencilRHSz, filterDecenteredStencilLHSz, filterDecenteredStencilRHSz, derivStartsZ, derivEndsZ, mesh.nz, [])

    # Add buffer zones into the derivative matrices
    if not hasattr(numMethods, 'changeOrderX') or numMethods.changeOrderX:
        LHSx, RHSx = addBufferToMatrix(LHSx, LHSxb, RHSx, RHSxb, mesh.nx, mesh.x.buffer)
    if not hasattr(numMethods, 'changeOrderY') or numMethods.changeOrderY:
        LHSy, RHSy = addBufferToMatrix(LHSy, LHSyb, RHSy, RHSyb, mesh.ny, mesh.y.buffer)
    if not hasattr(numMethods, 'changeOrderZ') or numMethods.changeOrderZ:
        LHSz, RHSz = addBufferToMatrix(LHSz, LHSzb, RHSz, RHSzb, mesh.nz, mesh.z.buffer)

    # Transform due to mesh stretching

    # Chose the spatial derivative method for the metric
    if not hasattr(numMethods, 'metricMethod') or numMethods.metricMethod is None:  # If unspecified, use SL4
        numMethods.metricMethod = 'SL4'
    elif numMethods.metricMethod == '':  # If left blank, use same as spatial derivatives
        numMethods.metricMethod = numMethods.spatialDerivs

    if mesh.nz == 4:  # Use SL4 for z if only 4 nodes are being used, due to its smaller stencil
        numMethods.metricMethodz = 'SL4'
    else:
        numMethods.metricMethodz = numMethods.metricMethod

    LHSx = applyMetric(LHSx, mesh.X, mesh.x, domain.xf, numMethods.metricMethod)
    LHSy = applyMetric(LHSy, mesh.Y, mesh.y, domain.yf, numMethods.metricMethod)
    LHSz = applyMetric(LHSz, mesh.Z, mesh.z, domain.zf, numMethods.metricMethodz)

    # Remove ends of filters if needed
    if numMethods.filter_borders:
        if hasattr(numMethods, 'filter_borders_start_x') and not numMethods.filter_borders_start_x:
            for i in range(len(f_LHSx)):
                f_LHSx[i][:5, :] = 0
                f_RHSx[i][:5, :] = 0
                f_LHSx[i][:5, :5] = np.eye(5)
                f_RHSx[i][:5, :5] = np.eye(5)
        if hasattr(numMethods, 'filter_borders_end_x') and not numMethods.filter_borders_end_x:
            for i in range(len(f_LHSx)):
                f_LHSx[i][-4:, :] = 0
                f_RHSx[i][-4:, :] = 0
                f_LHSx[i][-4:, -4:] = np.eye(5)
                f_RHSx[i][-4:, -4:] = np.eye(5)
        if hasattr(numMethods, 'filter_borders_start_y') and not numMethods.filter_borders_start_y:
            for i in range(len(f_LHSy)):
                f_LHSy[i][:5, :] = 0
                f_RHSy[i][:5, :] = 0
                f_LHSy[i][:5, :5] = np.eye(5)
                f_RHSy[i][:5, :5] = np.eye(5)
        if hasattr(numMethods, 'filter_borders_end_y') and not numMethods.filter_borders_end_y:
            for i in range(len(f_LHSy)):
                f_LHSy[i][-4:, :] = 0
                f_RHSy[i][-4:, :] = 0
                f_LHSy[i][-4:, -4:] = np.eye(5)
                f_RHSy[i][-4:, -4:] = np.eye(5)
        if hasattr(numMethods, 'filter_borders_start_z') and not numMethods.filter_borders_start_z:
            for i in range(len(f_LHSz)):
                f_LHSz[i][:5, :] = 0
                f_RHSz[i][:5, :] = 0
                f_LHSz[i][:5, :5] = np.eye(5)
                f_RHSz[i][:5, :5] = np.eye(5)
        if hasattr(numMethods, 'filter_borders_end_z') and not numMethods.filter_borders_end_z:
            for i in range(len(f_LHSz)):
                f_LHSz[i][-4:, :] = 0
                f_RHSz[i][-4:, :] = 0
                f_LHSz[i][-4:, -4:] = np.eye(5)
                f_RHSz[i][-4:, -4:] = np.eye(5)

    # Save to output structure

    matrices = {}

    matrices.x.types = typeMapX
    matrices.y.types = typeMapY
    matrices.z.types = typeMapZ

    matrices.x.LHS = LHSx
    matrices.x.RHS = RHSx
    matrices.y.LHS = LHSy
    matrices.y.RHS = RHSy
    matrices.z.LHS = LHSz
    matrices.z.RHS = RHSz

    matrices.x.fLHS = f_LHSx
    matrices.x.fRHS = f_RHSx
    matrices.y.fLHS = f_LHSy
    matrices.y.fRHS = f_RHSy
    matrices.z.fLHS = f_LHSz
    matrices.z.fRHS = f_RHSz

def makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, derivStarts, derivEnds, n, bufferInfo=None):
    """
    This function creates the matrices for each direction
    """

    # Check if the mesh is large enough for the stencil
    if n != 1 and n < 2 * len(centeredStencilRHS) - 1:
        raise ValueError(f"Mesh is not large enough for one of the stencils. It has {n} nodes but the stencil needs at least {2 * len(centeredStencilRHS) - 1}.")

    nTypes = len(derivStarts)

    LHS = [None] * nTypes
    RHS = [None] * nTypes

    if n == 1:  # If single point in this direction, set derivative to zero and filter to one
        if centeredStencilRHS[0] == 0:
            for i in range(nTypes):
                LHS[i] = np.ones((n, n))
                RHS[i] = np.zeros((n, n))
        else:
            for i in range(nTypes):
                LHS[i] = np.ones((n, n))
                RHS[i] = np.ones((n, n))
    else:
        # Create the base stencils
        LHS_base = sp.diags(centeredStencilLHS[0] * np.ones(n), 0)
        RHS_base = sp.diags(centeredStencilRHS[0] * np.ones(n), 0)

        invert_stencil = -1 if centeredStencilRHS[0] == 0 else 1

        for i in range(1, len(centeredStencilLHS)):
            inds = np.arange(n)[:, np.newaxis] + np.mod(np.arange(n) - 1 + i, n)
            LHS_base[inds] = centeredStencilLHS[i]
            inds = np.arange(n)[:, np.newaxis] + np.mod(np.arange(n) - 1 - i, n)
            LHS_base[inds] = centeredStencilLHS[i]

        for i in range(1, len(centeredStencilRHS)):
            inds = np.arange(n)[:, np.newaxis] + np.mod(np.arange(n) - 1 + i, n)
            RHS_base[inds] = centeredStencilRHS[i]
            inds = np.arange(n)[:, np.newaxis] + np.mod(np.arange(n) - 1 - i, n)
            RHS_base[inds] = invert_stencil * centeredStencilRHS[i]

        # Check if this is a buffer zone that needs to be upwind
        # and add that to the list of starts and ends so that the decentered
        # stencil is used
        if bufferInfo is not None:
            if hasattr(bufferInfo.i, 'upwind') and bufferInfo.i.upwind:
                for i in range(nTypes):
                    derivStarts[i] = np.unique(np.concatenate([derivStarts[i], np.arange(1, bufferInfo.i.n + 1)]), axis=0)[::-1]

            if hasattr(bufferInfo.f, 'upwind') and bufferInfo.f.upwind:
                for i in range(nTypes):
                    derivEnds[i] = np.unique(np.concatenate([derivEnds[i], np.arange(n - bufferInfo.f.n + 1, n + 1)]), axis=0)

        # Add startings and endings

        mLHS, nLHS = decenteredStencilLHS.shape
        mRHS, nRHS = decenteredStencilRHS.shape

        for i in range(nTypes):
            LHS_temp = LHS_base.copy()
            RHS_temp = RHS_base.copy()

            for j in range(len(derivStarts[i])):
                ind_start = derivStarts[i][j]

                LHS_temp[ind_start:ind_start + mLHS - 1, :] = 0
                RHS_temp[ind_start:ind_start + mRHS - 1, :] = 0

                LHS_temp[ind_start:ind_start + mLHS - 1, ind_start:ind_start + nLHS - 1] = decenteredStencilLHS
                RHS_temp[ind_start:ind_start + mRHS - 1, ind_start:ind_start + nRHS - 1] = decenteredStencilRHS

            for j in range(len(derivEnds[i])):
                ind_end = derivEnds[i][j]

                LHS_temp[ind_end - mLHS + 1:ind_end, :] = 0
                RHS_temp[ind_end - mRHS + 1:ind_end, :] = 0

                LHS_temp[ind_end - mLHS + 1:ind_end, ind_end - nLHS + 1:ind_end] = np.flip(decenteredStencilLHS, axis=0)
                RHS_temp[ind_end - mRHS + 1:ind_end, ind_end - nRHS + 1:ind_end] = invert_stencil * np.flip(decenteredStencilRHS, axis=0)

            LHS[i] = LHS_temp
            RHS[i] = RHS_temp

    return LHS, RHS

def addBufferToMatrix(baseMatrixL, bufferMatrixL, baseMatrixR, bufferMatrixR, n, bufferInfo):
    """
    This function adds the buffer zone to the matrices
    """

    ni = bufferInfo.i.n
    nf = bufferInfo.f.n

    if hasattr(bufferInfo.i, 'transition'):
        nti = int(bufferInfo.i.transition * ni)
    else:
        nti = ni
    if hasattr(bufferInfo.f, 'transition'):
        ntf = int(bufferInfo.f.transition * nf)
    else:
        ntf = nf

    ni1 = ni - nti + 1
    ni2 = ni
    nf1 = n - nf + 1
    nf2 = nf1 + ntf - 1

    # The buffer zone can be computed by a different type of derivatives. The transition is done smoothly.

    eta = np.ones((n, 1))

    if ni > 0:
        eta[:ni1] = 0
        eta[ni1:ni2] = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, ni2 - ni1 + 1))

    if nf > 0:
        eta[nf2:] = 0
        eta[nf1:nf2] = 0.5 + 0.5 * np.cos(np.linspace(0, np.pi, nf2 - nf1 + 1))

    nTypes = len(baseMatrixL)
    newMatrixL = [None] * nTypes
    newMatrixR = [None] * nTypes

    eta = np.sqrt(eta)

    for i in range(nTypes):
        newMatrixL[i] = eta * baseMatrixL[i] + (1 - eta) * bufferMatrixL[i]
        newMatrixR[i] = eta * baseMatrixR[i] + (1 - eta) * bufferMatrixR[i]

    return newMatrixL, newMatrixR

def applyMetric(LHS, X, meshInfo, xf, method):
    """
    This function applies the metric to the LHS matrix
    """

    # Scale the LHS matrix due to the mesh stretching. SL4 method is used here

    if meshInfo.periodic and meshInfo.fix_periodic_domain_size and len(X) > 1:  # Add temp node for a periodic mesh
        X = np.concatenate([X, [xf]])
        added_temp_node = True
    else:
        added_temp_node = False

    n = len(X)
    if n == 1:
        return

    centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS = finiteDifferenceCoefficients(method)
    LHS_temp, RHS_temp = makeMatricesEachDirection(centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS, [1], [n], n, [])

    dXdEta = LHS_temp @ X[:, np.newaxis]

    if added_temp_node:
        dXdEta[-1] = 0

    nTypes = len(LHS)
    for i in range(nTypes):
        LHS[i] = LHS[i] @ dXdEta


