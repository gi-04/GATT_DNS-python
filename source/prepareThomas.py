import numpy as np

def prepareThomas(matrix):
    """
    This function takes the matrices as inputs and prepares the vectors for the Fortran solver.
    Some computations that would be repeated during runtime have their results precomputed here.
    """
    
    # Get sizes
    isPeriodic = matrix['LHS'][0][-1, 0] != 0
    nTypes = np.max(matrix['types'])
    N = matrix['LHS'][0].shape[0]

    # Prepare LHS
    A = np.zeros((N-1, nTypes))
    B = np.zeros((N, nTypes))
    C = np.zeros((N-1, nTypes))
    D = np.zeros((N, nTypes))

    A1 = np.zeros((1, nTypes))
    Cn = np.zeros((1, nTypes))

    Af = np.zeros((N-1, nTypes))
    Bf = np.zeros((N, nTypes))
    Cf = np.zeros((N-1, nTypes))
    Df = np.zeros((N, nTypes))

    A1f = np.zeros((1, nTypes))
    Cnf = np.zeros((1, nTypes))

    for i in range(nTypes):
        A[:, i] = np.diag(matrix['LHS'][i], -1)
        B[:, i] = np.diag(matrix['LHS'][i])
        C[:, i] = np.diag(matrix['LHS'][i], 1)

        Af[:, i] = np.diag(matrix['fLHS'][i], -1)
        Bf[:, i] = np.diag(matrix['fLHS'][i])
        Cf[:, i] = np.diag(matrix['fLHS'][i], 1)

    if isPeriodic:
        for i in range(nTypes):
            A1[0, i] = matrix['LHS'][i][0, -1]
            Cn[0, i] = matrix['LHS'][i][-1, 0]

            A1f[0, i] = matrix['fLHS'][i][0, -1]
            Cnf[0, i] = matrix['fLHS'][i][-1, 0]

        A = A[1:, :]
        B = B[1:, :]
        C = C[1:, :]
        Af = Af[1:, :]
        Bf = Bf[1:, :]
        Cf = Cf[1:, :]

        for i in range(N-2):
            C[i, :] = C[i, :] / B[i, :]
            B[i+1, :] = B[i+1, :] - A[i, :] * C[i, :]

            Cf[i, :] = Cf[i, :] / Bf[i, :]
            Bf[i+1, :] = Bf[i+1, :] - Af[i, :] * Cf[i, :]

        B = 1.0 / B
        Bf = 1.0 / Bf

        # A1 and Cn will be stored as the first elements of A and C
        A = np.vstack([A1, A])
        Af = np.vstack([A1f, Af])
        C = np.vstack([Cn, C])
        Cf = np.vstack([Cnf, Cf])

        for i in range(nTypes):
            Dtemp = np.linalg.inv(matrix['LHS'][i])
            D[:, i] = Dtemp[0, :]

            Dtemp = np.linalg.inv(matrix['fLHS'][i])
            Df[:, i] = Dtemp[0, :]

    else:
        A1[:] = 0
        Cn[:] = 0
        A1f[:] = 0
        Cnf[:] = 0

        for i in range(N-1):
            C[i, :] = C[i, :] / B[i, :]
            B[i+1, :] = B[i+1, :] - A[i, :] * C[i, :]

            Cf[i, :] = Cf[i, :] / Bf[i, :]
            Bf[i+1, :] = Bf[i+1, :] - Af[i, :] * Cf[i, :]

        B = 1.0 / B
        Bf = 1.0 / Bf

    # Prepare RHS
    done = False

    RHSTemp = np.zeros((N, N, nTypes))

    for i in range(nTypes):
        RHSTemp[:, :, i] = matrix['RHS'][i]

    RHSDiag = fullDiag(RHSTemp, 0)

    nDiags = 1
    while not done:
        nextDiag = fullDiag(RHSTemp, nDiags)
        prevDiag = fullDiag(RHSTemp, -nDiags)

        if (np.any(nextDiag) or np.any(prevDiag)) and not (2 * nDiags - 1 > N):
            nDiags += 1
            RHSDiag = np.hstack([prevDiag, RHSDiag, nextDiag])
        else:
            done = True

    done = False

    RHSTemp = np.zeros((N, N, nTypes))

    for i in range(nTypes):
        RHSTemp[:, :, i] = matrix['fRHS'][i]

    RHSDiagf = fullDiag(RHSTemp, 0)

    nDiagsf = 1
    while not done:
        nextDiag = fullDiag(RHSTemp, nDiagsf)
        prevDiag = fullDiag(RHSTemp, -nDiagsf)

        if (np.any(nextDiag) or np.any(prevDiag)) and not (2 * nDiagsf - 1 > N):
            nDiagsf += 1
            RHSDiagf = np.hstack([prevDiag, RHSDiagf, nextDiag])
        else:
            done = True

    # Store outputs
    matrix['periodic'] = isPeriodic
    matrix['A'] = A
    matrix['B'] = B
    matrix['C'] = C
    matrix['D'] = D
    matrix['Af'] = Af
    matrix['Bf'] = Bf
    matrix['Cf'] = Cf
    matrix['Df'] = Df
    matrix['R'] = RHSDiag
    matrix['nRHS'] = nDiags
    matrix['Rf'] = RHSDiagf
    matrix['nRHSf'] = nDiagsf
    matrix['nTypes'] = nTypes

    return matrix

def fullDiag(M, k):
    """
    Extracts the k-th diagonal from each 2D slice of a 3D array M.
    """
    D = np.zeros((M.shape[0], 1, M.shape[2]))

    for i in range(M.shape[2]):
        D[:, 0, i] = np.diag(np.roll(M[:, :, i], -k, axis=1))

    return D



