import os
import numpy as np
import scipy.io
import subprocess
import re

def runBatch(parName, stage=None):
    if stage is None:
        stage = [1, 2, 3]

    if isinstance(parName, list):
        for i in range(len(parName)):
            runBatch(parName[i], stage)
        return

    parName = globals()[parName]
    maxNumCompThreads(p_row * p_col)

    if any(s == 1 for s in stage):
        info = None
        if 'info' in globals():
            info = globals()['info']
        else:
            mesh = scipy.io.loadmat(f"{caseName}/mesh.mat")['mesh']
        Ldomain = mesh.X[-1] - mesh.X[0]
        c = 1 / Ma
        v1 = c + 1
        v2 = c - 1
        T1 = Ldomain / v1
        T2 = Ldomain / v2
        Twave = 4 * (T1 + T2)
        eta = 0.7
        filesForMean = round(eta * Twave / time.qtimes)
        spacingForMean = round(1.1 * filesForMean)
        nMeans = 10
        tol = 1e-12
        runDNSautomean_internal(parName, caseName, time.qtimes, time.control, filesForMean, spacingForMean, nMeans, tol)

    if any(s == 2 for s in stage):
        mesh = scipy.io.loadmat(f"{caseName}/mesh.mat")['mesh']
        Ldomain = mesh.X[-1] - mesh.X[0]
        c = 1 / Ma
        v1 = c + 1
        v2 = c - 1
        T1 = Ldomain / v1
        T2 = Ldomain / v2
        Twave = 4 * (T1 + T2)
        eta = 0.7
        filesForMean = round(eta * Twave / time.qtimes)
        spacingForMean = round(1.1 * filesForMean)
        nMeans = 10
        tol = 1e-12
        runDNSautomean_internal(parName, caseName, time.qtimes, time.control, filesForMean, spacingForMean, nMeans, tol)

    if any(s == 3 for s in stage):
        calcInstability_internal(caseName)

def calcInstability_internal(caseFolder):
    # Define parameters
    comment = 'SFD01'

    # Method parameters
    M = 2000  # Size of Krylov space
    epsilon = 1e-6  # Disturbance magnitude
    dT = 2
    singleBeta = True
    firstOrder = True

    # SFD parameters (optional)
    resumeSFD = True

    # Domain parameters (optional)
    # nz = 16
    # beta = np.pi

    # Domain decomposition (optional)
    # global_p_row = 4
    # global_p_col = 1

    # Output parameters
    nModes = 50  # Number of modes to be saved
    removeCC = True
    nKrilov = M  # When the results files will be generated
    saveEvery = 100  # When intermediary files will be saved for resuming

    printEigsEvery = 10  # When intermediary eigenvalues will be printed
    printEigsN = 6  # Number of eigenvalues printed

    displayAllOutputs = False

    # Set folders for libraries
    matlabDir = ''  # Leave empty for automatic directory
    decompDir = '/usr/local/2decomp_fft'

    # Parameters required by the DNS
    logAll = False
    displayCompiling = False
    optimizeCode = True
    debugger = False
    profiler = False
    runningLST = True

    if displayAllOutputs:
        suppressOutput = ''
    else:
        suppressOutput = ' > /dev/null 2>&1'

    # Read parameters and change some values
    caseFile = 'parameters'
    if os.path.exists(f"{caseFolder}/{caseFile}.m"):
        with open(f"{caseFolder}/{caseFile}.m", "r") as f:
            code = f.read()
        exec(code)
    elif os.path.exists(f"{caseFolder}/bin/{caseFile}.m"):
        with open(f"{caseFolder}/bin/{caseFile}.m", "r") as f:
            code = f.read()
        exec(code)
    else:
        raise FileNotFoundError(f"Parameters file not found in {caseFolder}")

    caseName = caseFolder
    logAll = False

    if 'beta' in globals() and globals()['beta'] > 0:
        domain.zi = 0
        domain.zf = 2 * np.pi / globals()['beta']
        mesh.z.fixPeriodicDomainSize = True
        if singleBeta:
            mesh.z.n = 4
    if 'nz' in globals():
        mesh.z.n = globals()['nz']

    if 'SFD' in globals():
        extraVars = globals()['SFD'].keys()
        for var in extraVars:
            exec(f"numMethods.SFD.{var} = globals()['SFD']['{var}']")

    if 'global_p_row' in globals():
        p_row = globals()['global_p_row']
    if 'global_p_col' in globals():
        p_col = globals()['global_p_col']

    time.control = 'cfl'
    time.qtimes = dT
    time.tmax = dT

    # Create case name
    if mesh.z.n == 1:
        singleBeta = False
    caseNameInstability = f"{caseFolder}-dT{dT}-eps{epsilon}"
    if not firstOrder:
        caseNameInstability += '-SO'
    if mesh.z.n > 1 and not singleBeta:
        caseNameInstability += '-MB'
    if 'nz' in globals():
        caseNameInstability += f"-nz{globals()['nz']}"
    if 'beta' in globals():
        caseNameInstability += f"-beta{globals()['beta']/np.pi}pi"
    if 'comment' in globals() and globals()['comment'] != '':
        caseNameInstability += f"-{globals()['comment']}"
    print(f"Case name: {caseNameInstability}")

    # Run preprocessing
    print("Running preprocessor")
    import source
    import source.boundaries

    # Clear previous system calls cache
    subprocess.run([""], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Create a copy of the original bin folder and mesh file to be restored after compling
    cleanupObj = lambda: cleanUpFunction(caseFolder)
    subprocess.run(f"cp -r {caseFolder}/bin {caseFolder}/bin_backup {suppressOutput}", shell=True)
    subprocess.run(f"cp {caseFolder}/mesh.mat {caseFolder}/mesh.mat.backup {suppressOutput}", shell=True)

    preprocessing()

    # Set maximum number of threads for Matlab processes
    maxNumCompThreads(p_row * p_col)

    # Fix nSave in parameters.F90 and set to 0
    subprocess.run(f"sed -i 's/    integer :: nSave = .*/    integer :: nSave = 0/' {caseFolder}/bin/parameters.F90 {suppressOutput}", shell=True)

    # Fix resumeSFD in parameters.F90
    if resumeSFD:
        subprocess.run(f"sed -i 's/    integer :: resumeMeanFlow = .*/    integer :: resumeMeanFlow = 1/' {caseFolder}/bin/parameters.F90 {suppressOutput}", shell=True)

    # Compile code
    print("Compiling code")
    compileFortran()

    # Create new folders and copy files
    print("Creating new folders and files")

    if not os.path.exists(f"{caseFolder}/Instability"):
        instFolderName = f"/dev/shm/Instability{np.random.randint(1e8)}"
        os.makedirs(instFolderName)
        subprocess.run(f"ln -s {instFolderName} {caseFolder}/Instability {suppressOutput}", shell=True)
    if not os.path.exists(f"{caseFolder}/Instability/{caseNameInstability}"):
        os.makedirs(f"{caseFolder}/Instability/{caseNameInstability}")
    else:
        subprocess.run(f"rm {caseFolder}/Instability/{caseNameInstability}/* -r {suppressOutput}", shell=True)

    subprocess.run(f"cp {caseFolder}/bin {caseFolder}/Instability/{caseNameInstability}/bin -r {suppressOutput}", shell=True)

    logFile = open(f"{caseFolder}/Instability/{caseNameInstability}/log.txt", "w")
    logFile.close()

    # Restore the original files
    cleanupObj()
    clear_function(cleanupObj)

    # Create new cleanUpObj to remove the Instability folder
    cleanupObj = lambda: cleanUpFunction(caseFolder, caseNameInstability)

    # Read the base flow
    print("Reading base flow: ")

    if os.path.exists(f"{caseFolder}/baseflow.mat"):
        baseflow = scipy.io.loadmat(f"{caseFolder}/baseflow.mat")
        print(f"{caseFolder}/baseflow.mat")
    else:
        nStep = checkPreviousRun(caseFolder)
        baseflow = scipy.io.loadmat(f"{caseFolder}/flow_{nStep:010d}.mat")
        print(f"{caseFolder}/flow_{nStep:010d}.mat")

    # If needed, replicate flow to make it 3D
    if mesh.z.n > 1 and baseflow['U'].shape[2] == 1:
        baseflow['U'] = np.repeat(baseflow['U'][:, :, np.newaxis], mesh.z.n, axis=2)
        baseflow['V'] = np.repeat(baseflow['V'][:, :, np.newaxis], mesh.z.n, axis=2)
        baseflow['W'] = np.repeat(baseflow['W'][:, :, np.newaxis], mesh.z.n, axis=2)
        baseflow['R'] = np.repeat(baseflow['R'][:, :, np.newaxis], mesh.z.n, axis=2)
        baseflow['E'] = np.repeat(baseflow['E'][:, :, np.newaxis], mesh.z.n, axis=2)

    # Get the physical flow region
    flowRegion = ~np.isnan(baseflow['U'])  # Get region outside of walls

    if removeBufferZone:
        flowRegion[[0:mesh.x.buffer.i.n, -mesh.x.buffer.f.n:-1], :, :] = False  # Remove buffer zones
        flowRegion[:, [0:mesh.y.buffer.i.n, -mesh.y.buffer.f.n:-1], :] = False
        flowRegion[:, :, [0:mesh.z.buffer.i.n, -mesh.z.buffer.f.n:-1]] = False

    if mesh.z.n == 1:
        N = 4 * np.sum(flowRegion)
    elif not singleBeta:
        N = 5 * np.sum(flowRegion)
    else:
        N = 5 * np.sum(np.sum(flowRegion[:, :, 0], axis=0))

    # Allocate variables
    print("Allocating variables")
    zeta = np.zeros((N, M))
    H = np.zeros((M, M))

    # Check for previous results and reload data
    if os.path.exists(f"{caseFolder}/results-{caseNameInstability}.mat"):
        print("Previous results file found, resuming")
        previousResults = scipy.io.loadmat(f"{caseFolder}/results-{caseNameInstability}.mat")

        M0 = min(previousResults['M'], M)
        zeta[:, 0:M0] = previousResults['zeta'][:, 0:M0]
        H[0:M0, 0:M0 - 1] = previousResults['H'][0:M0, 0:M0 - 1]
    else:
        M0 = 1

        zeta[:, 0] = calcInitialDisturb(mesh.X, mesh.Y, mesh.Z, singleBeta, flowRegion)

    # Set the initial SFD state
    if resumeSFD:
        t = 0
        U = baseflow['U']
        V = baseflow['V']
        W = baseflow['W']
        R = baseflow['R']
        E = baseflow['E']
        scipy.io.savemat(f"{caseFolder}/Instability/{caseNameInstability}/meanflowSFD.mat", {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t}, oned_as='column', do_compression=True)

    # Iterate
    if firstOrder:  # If using a first order approximation, compute the baseflow drift
        print("Iteration 0")
        t = 0

        U = baseflow['U']
        V = baseflow['V']
        W = baseflow['W']
        R = baseflow['R']
        E = baseflow['E']
        scipy.io.savemat(f"{caseFolder}/Instability/{caseNameInstability}/flow_0000000000.mat", {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t}, oned_as='column', do_compression=True)

        subprocess.run(f"(cd {caseFolder}/Instability/{caseNameInstability}/bin && mpirun -np {p_row*p_col} main {caseNameInstability}) {suppressOutput}", shell=True)
        flowMinus = scipy.io.loadmat(f"{caseFolder}/Instability/{caseNameInstability}/flow_0000000001.mat")
        os.remove(f"{caseFolder}/Instability/{caseNameInstability}/flow_0000000001.mat")

        UM = flow2vec(flowMinus, flowRegion, singleBeta)

    tic = time.time()
    elapsedTime = []
    for k in range(M0, M):
        print("\n")
        print(f"Iteration {k} of {M}")

        # Get current disturbance
        disturbance = vec2flow(zeta[:, k] * epsilon * np.sqrt(N), flowRegion, singleBeta)

        # Save flows, run DNS and read results
        t = 0

        U = baseflow['U'] + disturbance['U']
        V = baseflow['V'] + disturbance['V']
        W = baseflow['W'] + disturbance['W']
        R = baseflow['R'] + disturbance['R']
        E = baseflow['E'] + disturbance['E']
        scipy.io.savemat(f"{caseFolder}/Instability/{caseNameInstability}/flow_0000000000.mat", {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t}, oned_as='column', do_compression=True)

        subprocess.run(f"(cd {caseFolder}/Instability/{caseNameInstability}/bin && mpirun -np {p_row*p_col} main {caseNameInstability}) {suppressOutput}", shell=True)
        flowPlus = scipy.io.loadmat(f"{caseFolder}/Instability/{caseNameInstability}/flow_0000000001.mat")
        os.remove(f"{caseFolder}/Instability/{caseNameInstability}/flow_0000000001.mat")

        UP = flow2vec(flowPlus, flowRegion, singleBeta)

        if not firstOrder:
            U = baseflow['U'] - disturbance['U']
            V = baseflow['V'] - disturbance['V']
            W = baseflow['W'] - disturbance['W']
            R = baseflow['R'] - disturbance['R']
            E = baseflow['E'] - disturbance['E']
            scipy.io.savemat(f"{caseFolder}/Instability/{caseNameInstability}/flow_0000000000.mat", {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t}, oned_as='column', do_compression=True)

            subprocess.run(f"(cd {caseFolder}/Instability/{caseNameInstability}/bin && mpirun -np {p_row*p_col} main {caseNameInstability}) {suppressOutput}", shell=True)
            flowMinus = scipy.io.loadmat(f"{caseFolder}/Instability/{caseNameInstability}/flow_0000000001.mat")
            os.remove(f"{caseFolder}/Instability/{caseNameInstability}/flow_0000000001.mat")

            UM = flow2vec(flowMinus, flowRegion, singleBeta)

            # Compute matrices
            B = (UP - UM) / (2 * epsilon * np.sqrt(N))
        else:
            B = (UP - UM) / (epsilon * np.sqrt(N))

        uPrime = B

        for j in range(1, k + 1):
            H[j, k] = np.dot(zeta[:, j], B)
            uPrime -= H[j, k] * zeta[:, j]

        if k < M - 1:
            H[k + 1, k] = np.linalg.norm(uPrime)
            zeta[:, k + 1] = uPrime / H[k + 1, k]

        if np.any(np.isnan(H)) or np.any(np.isinf(H)):
            raise ValueError("DNS solution has failed")

        # Save results to file
        if (k + 1) % saveEvery == 0 and k != M0:
            actualM = M
            actualH = H
            actualZeta = zeta

            M = k + 1

            H[M:, :] = []
            H[:, M:] = []
            zeta[:, M:] = []

            print(f"Saving results to results-{caseNameInstability}.mat")
            scipy.io.savemat(f"{caseFolder}/results-{caseNameInstability}.mat", {'zeta': zeta, 'H': H, 'M': M}, oned_as='column', do_compression=True)

            zeta = actualZeta
            H = actualH
            M = actualM

        # Compute modes
        if (k + 1) in nKrilov:
            print("Computing modes")

            psiH, lambdaH = np.linalg.eig(H[:k + 1, :k + 1])
            lambdaB = np.diag(lambdaH)

            lambdaA = np.log(lambdaB) / dT
            index = np.argsort(np.real(lambdaA))[::-1]
            lambdaA = lambdaA[index]

            psiH = psiH[:, index]

            psiA = np.dot(zeta[:, :k + 1], psiH[:, :min(nModes, k + 1)])

            modes = vec2modes(psiA, flowRegion, singleBeta, min(nModes, k + 1), lambdaA, removeCC)
            modes['lambda'] = lambdaA
            modes['X'] = mesh.X
            modes['Y'] = mesh.Y
            modes['Z'] = mesh.Z

            # Save modes
            print(f"Saving modes to modes-{caseNameInstability}-M{k}.mat")
            scipy.io.savemat(f"{caseFolder}/modes-{caseNameInstability}-M{k}.mat", {'modes': modes}, oned_as='column', do_compression=True)

        # Print remaining time
        elapsedTime.append(time.time() - tic)
        remainingTime = (np.mean(elapsedTime[1:]) * (M - k - 1))
        if len(elapsedTime) > 100:
            elapsedTime.pop(0)
        print(f"Elapsed time: {sec2str(elapsedTime[-1])}")
        print(f"Estimated time remaining: {sec2str(remainingTime)}")

        # Print leading eigenvalues
        if (k + 1) % printEigsEvery == 0:
            leadingEigs = scipy.sparse.linalg.eigs(H[:k + 1, :k + 1], k=min(k + 1, printEigsN), which='LM', tol=1e-4, maxiter=max(300, 2 * k + 2))
            leadingEigs = np.log(leadingEigs) / dT
            print(f"Leading eigenvalues : {leadingEigs}")
            if 'leadingEigsPrev' in globals():
                print(f"Change : {np.abs((leadingEigs[:len(leadingEigsPrev)] - leadingEigsPrev) / leadingEigs[:len(leadingEigsPrev)])}")
            leadingEigsPrev = leadingEigs

def cleanUpFunction(caseFolder, caseNameInstability=None):
    if caseNameInstability is None:
        subprocess.run(f"rm {caseFolder}/Instability {suppressOutput}", shell=True)
    else:
        subprocess.run(f"rm -r {caseFolder}/Instability/{caseNameInstability} {suppressOutput}", shell=True)

def preprocessing():
    # Preprocessing code goes here
    pass

def compileFortran():
    # Compile Fortran code
    subprocess.run(f"gfortran -o main *.f {suppressOutput}", shell=True)

def checkPreviousRun(caseFolder):
    # Check for previous run and return the latest step number
    flowFiles = [f for f in os.listdir(caseFolder) if f.startswith('flow_') and f.endswith('.mat')]
    flowFiles.sort(key=lambda x: int(x[5:-4]))
    return int(flowFiles[-1][5:-4])

def flow2vec(flow, flowRegion, singleBeta):
    # Convert flow to vector
    U = flow['U'][flowRegion]
    V = flow['V'][flowRegion]
    W = flow['W'][flowRegion]
    R = flow['R'][flowRegion]
    E = flow['E'][flowRegion]
    if singleBeta:
        return np.concatenate((U.flatten(), V.flatten(), R.flatten(), E.flatten()))
    else:
        return np.concatenate((U.flatten(), V.flatten(), W.flatten(), R.flatten(), E.flatten()))

def vec2flow(zeta, flowRegion, singleBeta):
    # Convert vector to flow
    N = zeta.size
    if singleBeta:
        U = np.reshape(zeta[:N // 4], flowRegion[0].shape)
        V = np.reshape(zeta[N // 4:N // 2], flowRegion[0].shape)
        R = np.reshape(zeta[N // 2:3 * N // 4], flowRegion[0].shape)
        E = np.reshape(zeta[3 * N // 4:], flowRegion[0].shape)
        flow = {'U': U, 'V': V, 'W': np.zeros_like(U), 'R': R, 'E': E}
    else:
        U = np.reshape(zeta[:N // 5], flowRegion[0].shape)
        V = np.reshape(zeta[N // 5:2 * N // 5], flowRegion[0].shape)
        W = np.reshape(zeta[2 * N // 5:3 * N // 5], flowRegion[0].shape)
        R = np.reshape(zeta[3 * N // 5:4 * N // 5], flowRegion[0].shape)
        E = np.reshape(zeta[4 * N // 5:], flowRegion[0].shape)
        flow = {'U': U, 'V': V, 'W': W, 'R': R, 'E': E}
    flow['U'][~flowRegion] = np.nan
    flow['V'][~flowRegion] = np.nan
    flow['W'][~flowRegion] = np.nan
    flow['R'][~flowRegion] = np.nan
    flow['E'][~flowRegion] = np.nan
    return flow

def vec2modes(psiA, flowRegion, singleBeta, nModes, lambdaA, removeCC):
    # Convert vector to modes
    modes = {'psi': np.zeros((nModes,) + flowRegion[0].shape, dtype=np.complex128)}
    for i in range(nModes):
        modes['psi'][i] = psiA[:, i].reshape(flowRegion[0].shape)
        if removeCC:
            modes['psi'][i] -= np.mean(modes['psi'][i])
    modes['lambda'] = lambdaA
    modes['X'] = mesh.X
    modes['Y'] = mesh.Y
    modes['Z'] = mesh.Z
    return modes

def sec2str(sec):
    # Convert seconds to formatted string
    h = int(sec // 3600)
    m = int((sec % 3600) // 60)
    s = int(sec % 60)
    return f"{h:02d}:{m:02d}:{s:02d}"

def calcInitialDisturb(X, Y, Z, singleBeta, flowRegion):
    # Calculate initial disturbance
    if singleBeta:
        disturbance = np.random.randn(*flowRegion[0].shape) + 1j * np.random.randn(*flowRegion[0].shape)
    else:
        disturbance = np.random.randn(*flowRegion[0].shape) + 1j * np.random.randn(*flowRegion[0].shape)
        disturbance = np.repeat(disturbance[np.newaxis, :, :], mesh.z.n, axis=0)
    disturbance /= np.linalg.norm(disturbance)
    return disturbance.flatten()

def flow2vec(flow, flowRegion, singleBeta):
    threeD = len(flow.shape) > 2
    
    if not threeD:  # 2D flow
        vec = np.concatenate([flow['U'][flowRegion], flow['V'][flowRegion], flow['R'][flowRegion], flow['E'][flowRegion]])
        
    elif not singleBeta:  # Fully 3D flow
        vec = np.concatenate([flow['U'][flowRegion], flow['V'][flowRegion], flow['W'][flowRegion], flow['R'][flowRegion], flow['E'][flowRegion]])
        
    else:  # Single beta, 3D flow
        if flow['U'].shape[2] > 1:
            nz = flow['U'].shape[2]
            Ztemp = np.linspace(0, 2 * np.pi, nz + 1)
            Ztemp = np.delete(Ztemp, -1)
            Ztemp = np.moveaxis(Ztemp, 0, -1)

            cosZ = np.cos(Ztemp) * 2 / nz
            sinZ = np.sin(Ztemp) * 2 / nz
            
            flow['U'] = np.sum(cosZ * flow['U'], axis=2)
            flow['V'] = np.sum(cosZ * flow['V'], axis=2)
            flow['W'] = np.sum(sinZ * flow['W'], axis=2)
            flow['R'] = np.sum(cosZ * flow['R'], axis=2)
            flow['E'] = np.sum(cosZ * flow['E'], axis=2)
        
        flowRegion = flowRegion[:, :, 0]
        
        vec = np.concatenate([flow['U'][flowRegion], flow['V'][flowRegion], flow['W'][flowRegion], flow['R'][flowRegion], flow['E'][flowRegion]])
    
    return vec

def vec2flow(vec, flowRegion, singleBeta):
    [nx, ny, nz] = flowRegion.shape
    N = np.sum(flowRegion)
    
    threeD = nz > 1
    
    if not threeD:  # 2D flow
        flow = {'U': np.zeros((nx, ny)), 'V': np.zeros((nx, ny)), 'W': np.zeros((nx, ny)), 'R': np.zeros((nx, ny)), 'E': np.zeros((nx, ny))}
        
        flow['U'][flowRegion] = vec[:N]
        flow['V'][flowRegion] = vec[N:2 * N]
        flow['R'][flowRegion] = vec[2 * N:3 * N]
        flow['E'][flowRegion] = vec[3 * N:4 * N]
        
    elif not singleBeta:  # Fully 3D flow
        flow = {'U': np.zeros((nx, ny, nz)), 'V': np.zeros((nx, ny, nz)), 'W': np.zeros((nx, ny, nz)), 'R': np.zeros((nx, ny, nz)), 'E': np.zeros((nx, ny, nz))}
        
        flow['U'][flowRegion] = vec[:N]
        flow['V'][flowRegion] = vec[N:2 * N]
        flow['W'][flowRegion] = vec[2 * N:3 * N]
        flow['R'][flowRegion] = vec[3 * N:4 * N]
        flow['E'][flowRegion] = vec[4 * N:5 * N]
    
    else:  # Single beta, 3D flow
        flow = {'U': np.zeros((nx, ny)), 'V': np.zeros((nx, ny)), 'W': np.zeros((nx, ny)), 'R': np.zeros((nx, ny)), 'E': np.zeros((nx, ny))}
        
        N = int(N / nz)
        
        flow['U'][flowRegion[:, :, 0]] = vec[:N]
        flow['V'][flowRegion[:, :, 0]] = vec[N:2 * N]
        flow['W'][flowRegion[:, :, 0]] = vec[2 * N:3 * N]
        flow['R'][flowRegion[:, :, 0]] = vec[3 * N:4 * N]
        flow['E'][flowRegion[:, :, 0]] = vec[4 * N:5 * N]
    
        Ztemp = np.linspace(0, 2 * np.pi, nz + 1)
        Ztemp = np.delete(Ztemp, -1)
        Ztemp = np.moveaxis(Ztemp, 0, -1)
        
        cosZ = np.cos(Ztemp)
        sinZ = np.sin(Ztemp)
        
        flow['U'] = flow['U'] * cosZ
        flow['V'] = flow['V'] * cosZ
        flow['W'] = flow['W'] * sinZ
        flow['R'] = flow['R'] * cosZ
        flow['E'] = flow['E'] * cosZ
    
    return flow

def vec2modes(vec, flowRegion, singleBeta, nModes, lambda_, removeCC):
    [nx, ny, nz] = flowRegion.shape
    N = np.sum(flowRegion)
    
    threeD = nz > 1
    
    if not threeD:  # 2D flow
        modes = {'U': np.full((nx, ny, 1, nModes), np.nan), 'V': np.full((nx, ny, 1, nModes), np.nan), 'W': np.full((nx, ny, 1, nModes), np.nan), 'R': np.full((nx, ny, 1, nModes), np.nan), 'E': np.full((nx, ny, 1, nModes), np.nan)}
        
        U = np.full((nx, ny), np.nan)
        V = np.full((nx, ny), np.nan)
        R = np.full((nx, ny), np.nan)
        E = np.full((nx, ny), np.nan)
        
        for i in range(nModes):
            if np.imag(lambda_[i]) >= 0:
                vec[:, i] = vec[:, i] / np.max(np.abs(vec[:, i]))
                
                U[flowRegion] = vec[:N, i]
                V[flowRegion] = vec[N:2 * N, i]
                R[flowRegion] = vec[2 * N:3 * N, i]
                E[flowRegion] = vec[3 * N:4 * N, i]
                
                modes['U'][:, :, 0, i] = U
                modes['V'][:, :, 0, i] = V
                modes['R'][:, :, 0, i] = R
                modes['E'][:, :, 0, i] = E
        
    elif not singleBeta:  # Fully 3D flow
        modes = {'U': np.full((nx, ny, nz, nModes), np.nan), 'V': np.full((nx, ny, nz, nModes), np.nan), 'W': np.full((nx, ny, nz, nModes), np.nan), 'R': np.full((nx, ny, nz, nModes), np.nan), 'E': np.full((nx, ny, nz, nModes), np.nan)}
        
        U = np.full((nx, ny, nz), np.nan)
        V = np.full((nx, ny, nz), np.nan)
        W = np.full((nx, ny, nz), np.nan)
        R = np.full((nx, ny, nz), np.nan)
        E = np.full((nx, ny, nz), np.nan)
        
        for i in range(nModes):
            if np.imag(lambda_[i]) >= 0:
                vec[:, i] = vec[:, i] / np.max(np.abs(vec[:, i]))
                
                U[flowRegion] = vec[:N, i]
                V[flowRegion] = vec[N:2 * N, i]
                W[flowRegion] = vec[2 * N:3 * N, i]
                R[flowRegion] = vec[3 * N:4 * N, i]
                E[flowRegion] = vec[4 * N:5 * N, i]
                
                modes['U'][:, :, :, i] = U
                modes['V'][:, :, :, i] = V
                modes['W'][:, :, :, i] = W
                modes['R'][:, :, :, i] = R
                modes['E'][:, :, :, i] = E
    
    else:  # Single beta, 3D flow
        N = int(N / nz)
        flowRegion = flowRegion[:, :, 0]
        
        modes = {'U': np.full((nx, ny, 1, nModes), np.nan), 'V': np.full((nx, ny, 1, nModes), np.nan), 'W': np.full((nx, ny, 1, nModes), np.nan), 'R': np.full((nx, ny, 1, nModes), np.nan), 'E': np.full((nx, ny, 1, nModes), np.nan)}
        
        U = np.full((nx, ny), np.nan)
        V = np.full((nx, ny), np.nan)
        W = np.full((nx, ny), np.nan)
        R = np.full((nx, ny), np.nan)
        E = np.full((nx, ny), np.nan)
        
        for i in range(nModes):
            if np.imag(lambda_[i]) >= 0:
                vec[:, i] = vec[:, i] / np.max(np.abs(vec[:, i]))
                
                U[flowRegion] = vec[:N, i]
                V[flowRegion] = vec[N:2 * N, i]
                W[flowRegion] = vec[2 * N:3 * N, i]
                R[flowRegion] = vec[3 * N:4 * N, i]
                E[flowRegion] = vec[4 * N:5 * N, i]
                
                modes['U'][:, :, 0, i] = U
                modes['V'][:, :, 0, i] = V
                modes['W'][:, :, 0, i] = W
                modes['R'][:, :, 0, i] = R
                modes['E'][:, :, 0, i] = E
    
    return modes

def calc_initial_disturb(X, Y, Z, singleBeta, flowRegion):
    x0_idx = np.argmin(np.gradient(X))
    y0_idx = np.argmin(np.gradient(Y))
    z0_idx = np.argmin(np.gradient(Z))
    
    x0 = X[x0_idx]
    y0 = Y[y0_idx]
    z0 = Z[z0_idx]
    
    alphaX = 10 / (X[-1] - X[0]) ** 2
    alphaY = 10 / (Y[-1] - Y[0]) ** 2
    alphaZ = 10 / (Z[-1] - Z[0]) ** 2
    
    if len(Z) == 1:
        eta = np.exp(-(alphaX * (X - x0) ** 2 + alphaY * (Y - y0) ** 2))
        
        flow = {'U': eta, 'V': eta, 'W': eta, 'R': eta, 'E': eta}
        
        zeta = flow2vec(flow, flowRegion, singleBeta)
        
    elif not singleBeta:
        eta = np.exp(-(alphaX * (X - x0) ** 2 + alphaY * (Y - y0) ** 2 + alphaZ * (np.moveaxis(Z, 0, -1) - z0) ** 2))
        
        flow = {'U': eta, 'V': eta, 'W': eta, 'R': eta, 'E': eta}
        
        zeta = flow2vec(flow, flowRegion, singleBeta)
    
    else:
        eta = np.exp(-(alphaX * (X - x0) ** 2 + alphaY * (Y - y0) ** 2))
        
        flow = {'U': eta, 'V': eta, 'W': eta, 'R': eta, 'E': eta}
        
        zeta = flow2vec(flow, flowRegion, singleBeta)
    
    zeta = zeta / np.linalg.norm(zeta)
    
    return zeta

def sec2str(sec):
    if sec < 60:
        return f"{sec}s"
    elif sec < 3600:
        min = sec // 60
        sec = sec % 60
        return f"{min}min {sec:02d}s"
    elif sec < 86400:
        hour = sec // 3600
        min = (sec % 3600) // 60
        sec = sec % 60
        return f"{hour}h {min:02d}min {sec:02d}s"
    else:
        day = sec // 86400
        hour = (sec % 86400) // 3600
        min = (sec % 3600) // 60
        sec = sec % 60
        return f"{day}d {hour:02d}h {min:02d}min {sec:02d}s"

def unpack_struct(structure):
    varList = []
    varNames = [varName[0] for varName in structure._fieldnames]
    
    for i, varName in enumerate(varNames):
        if isinstance(structure[varName], dict):
            varList2 = unpack_struct(structure[varName])
            for j, var in enumerate(varList2):
                varList.append(f"{varName}.{var}")
        else:
            varList.append(varName)
    
    return varList

def clean_up_function(caseFolder, caseNameInstability=None):
    if caseNameInstability is None:
        os.system(f"mv {caseFolder}/mesh.mat.backup {caseFolder}/mesh.mat >/dev/null 2>&1")
        os.system(f"rm -r {caseFolder}/bin >/dev/null 2>&1")
        os.system(f"mv {caseFolder}/bin_backup {caseFolder}/bin >/dev/null 2>&1")
    else:
        os.system(f"rm -r {caseFolder}/Instability/{caseNameInstability}")
        temp = os.listdir(f"{caseFolder}/Instability")
        if len(temp) == 2:
            os.system(f"rm -r $(readlink -f {caseFolder}/Instability)")
            os.system(f"rm -r {caseFolder}/Instability")

def run_dnsautomean_internal(parName, caseName, qtimes, timecontrol, filesForMean, spacingForMean, nMeans, tol):
    import re

    nMean = 1

    while True:
        changesStr = os.popen(f"tail -1 {caseName}/log.txt").read().strip()
        changes = [float(x) for x in re.findall(r'\d+\.?\d*', changesStr)]
        change = max(changes[5:10])
        
        if nMean > nMeans or change < tol:
            break
        
        print(f"Computing means loop: {nMean}\n Current change: {change}")
        nMean += 1

        caseFiles = []
        allFiles = os.listdir(caseName)
        for i, name in enumerate(allFiles):
            if len(name) == 19 and re.match(r'flow_\d*.mat', name):
                caseFiles.append(name)
        
        caseFiles = caseFiles[-filesForMean:]
        
        lastStep = int(caseFiles[-1][5:-4])
        
        U = np.zeros((1, 1))
        V = np.zeros((1, 1))
        W = np.zeros((1, 1))
        R = np.zeros((1, 1))
        E = np.zeros((1, 1))
        
        for i, caseFile in enumerate(caseFiles):
            print(caseFile)
            current = scipy.io.loadmat(f"{caseName}/{caseFile}")
            
            U += current['U']
            V += current['V']
            W += current['W']
            R += current['R']
            E += current['E']
        
        U /= filesForMean
        V /= filesForMean
        W /= filesForMean
        R /= filesForMean
        E /= filesForMean
        t = current['t'][0, 0]
        
        tstr = f"{caseName}/flow_{lastStep + 1:010d}.mat"
        scipy.io.savemat(tstr, {'t': t, 'U': U, 'V': V, 'W': W, 'R': R, 'E': E})
        
        os.system(f"mv {caseName}/meanflowSFD.mat {caseName}/meanflowSFD_{lastStep}.mat 2>/dev/null")
        
        if timecontrol == 'cfl':
            extraPars['time']['tmax'] = t + spacingForMean * qtimes
        elif timecontrol == 'dt':
            extraPars['time']['tmax'] = lastStep + spacingForMean * qtimes
        
        logFile2 = open(f"{caseName}/bin/log2.txt", "a")
        logFile2.write(f"Mean flow calculated from {caseFiles[0]} to {caseFiles[-1]}\n")
        logFile2.close()
        
        run_dns(parName, extraPars)


