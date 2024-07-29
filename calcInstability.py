import numpy as np
import os
import glob
import scipy.io

def calcInstability():
    for beta in np.arange(0, 1.1, 0.1) * np.pi:
        # Define parameters
        comment = ''
        caseFolder = '/home/felipe/autolst/teste'
        
        # Method parameters
        M = 8000  # Size of Krylov space
        epsilon = 1e-7  # Disturbance magnitude
        dT = 5
        singleBeta = True
        firstOrder = True
        removeBufferZone = False
        
        # SFD parameters (optional)
        resumeSFD = False
        SFD = {
            'type': 2,  # 0 = off, 1 = whole domain, 2 = buffer zone only;
            'Delta': float('inf'),
            'X': 0.005
        }
        
        # Domain parameters (optional)
        # nz = 16
        # beta = np.pi
        
        # Domain decomposition (optional)
        global_p_col = 4 if beta != 0 else 1
        p_row_max = 64 / global_p_col  # Number of cores/threads 
        global_p_row = get_p_row('parametersBaseFlow', p_row_max)

        # Output parameters
        nModes = [40, 60, 200, 400]  # Vector of modes to be saved
        removeCC = True
        nKrilov = [8000]  # When the results files will be generated
        saveEvery = 500  # When intermediary files will be saved for resuming 
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

        suppressOutput = '' if displayAllOutputs else ' > /dev/null 2>&1'

        # Read parameters and change some values
        caseFile = 'parameters'
        print('Reading parameters file: ')
        if os.path.isfile(os.path.join(caseFolder, f'{caseFile}.m')):  # First check case dir
            print(os.path.join(caseFolder, f'{caseFile}.m'))
            os.chdir(caseFolder)
            exec(open(caseFile).read())
            os.chdir('..')
        elif os.path.isfile(os.path.join(caseFolder, 'bin', f'{caseFile}.m')):  # Then check bin dir
            print(os.path.join(caseFolder, 'bin', f'{caseFile}.m'))
            os.chdir(caseFolder)
            os.chdir('bin')
            exec(open(caseFile).read())
            os.chdir('..')
            os.chdir('..')
        else:  # Now check for any .m in case dir
            allFiles = os.listdir(caseFolder)
            parFiles = [name for name in allFiles if name.endswith('.m')]
            
            if len(parFiles) == 1:
                print(os.path.join(caseFolder, parFiles[0]))
                os.chdir(caseFolder)
                exec(open(parFiles[0][:-2]).read())
                os.chdir('..')
            else:  # Finally, check for any .m in bin dir
                allFiles = os.listdir(os.path.join(caseFolder, 'bin'))
                parFiles = [name for name in allFiles if name.endswith('.m')]
                
                if len(parFiles) == 1:
                    print(os.path.join(caseFolder, 'bin', parFiles[0]))
                    os.chdir(caseFolder)
                    os.chdir('bin')
                    exec(open(parFiles[0][:-2]).read())
                    os.chdir('..')
                    os.chdir('..')
                else:
                    raise FileNotFoundError('Parameters file not found')

        caseName = caseFolder
        logAll = False

        if 'beta' in locals() and beta > 0:
            domain = {'zi': 0, 'zf': 2 * np.pi / beta}
            mesh = {'z': {'fixPeriodicDomainSize': True, 'n': 4 if singleBeta else None}}

        if 'nz' in locals():
            mesh['z']['n'] = nz

        if 'SFD' in locals():
            extraVars = unpackStruct(SFD)
            for var in extraVars:
                exec(f'numMethods.SFD.{var} = SFD.{var}')

        if 'global_p_row' in locals():
            p_row = global_p_row
        if 'global_p_col' in locals():
            p_col = global_p_col

        time = {'control': 'cfl', 'qtimes': dT, 'tmax': dT}

        # Create case name
        if mesh['z']['n'] == 1:
            singleBeta = False

        caseNameInstability = f"{caseFolder}-dT{dT}-eps{epsilon}"
        if not firstOrder:
            caseNameInstability += '-SO'
        if mesh['z']['n'] > 1 and not singleBeta:
            caseNameInstability += '-MB'
        if 'nz' in locals():
            caseNameInstability += f'-nz{nz}'
        if 'beta' in locals():
            caseNameInstability += f'-beta{beta/np.pi}pi'
        if 'comment' in locals() and comment:
            caseNameInstability += f'-{comment}'
        print(f'Case name: {caseNameInstability}')

        # Run preprocessing
        print('Running preprocessor')
        os.system('addpath source')
        os.system('addpath source/boundaries')

        # Set maximum number of threads for Python processes
        # maxNumCompThreads(p_row * p_col)  # Not applicable in Python

        # Clear previous system calls cache
        os.system('')

        # Create a copy of the original bin folder and mesh file to be restored after compiling
        cleanupObj = onCleanup(lambda: cleanUpFunction(caseFolder))
        os.system(f'cp -r {caseFolder}/bin {caseFolder}/bin_backup >/dev/null 2>&1')
        os.system(f'cp {caseFolder}/mesh.mat {caseFolder}/mesh.mat.backup >/dev/null 2>&1')

        preprocessing()

        # Fix nSave in parameters.F90 and set to 0
        os.system(f"sed -i 's/    integer :: nSave = .*/    integer :: nSave = 0/' {caseFolder}/bin/parameters.F90")

        # Fix resumeSFD in parameters.F90
        if resumeSFD:
            os.system(f"sed -i 's/    integer :: resumeMeanFlow = .*/    integer :: resumeMeanFlow = 1/' {caseFolder}/bin/parameters.F90")

        # Compile code
        print('Compiling code')
        compileFortran()

        # Create new folders and copy files
        print('Creating new folders and files')

        if not os.path.exists(os.path.join(caseFolder, 'Instability')):
            instFolderName = f"/dev/shm/Instability{np.random.randint(1e8)}"
            os.makedirs(instFolderName)
            os.symlink(instFolderName, os.path.join(caseFolder, 'Instability'))

        if not os.path.exists(os.path.join(caseFolder, 'Instability', caseNameInstability)):
            os.makedirs(os.path.join(caseFolder, 'Instability', caseNameInstability))
        else:
            os.system(f'rm {os.path.join(caseFolder, "Instability", caseNameInstability)}/* -r')

        os.system(f'cp {caseFolder}/bin {os.path.join(caseFolder, "Instability", caseNameInstability, "bin")} -r')

        logFile = open(os.path.join(caseFolder, 'Instability', caseNameInstability, 'log.txt'), 'w')
        logFile.close()

        # Restore the original files
        del cleanupObj
        del runningSFD

        # Create new cleanUpObj to remove the Instability folder
        cleanupObj = onCleanup(lambda: cleanUpFunction(caseFolder, caseNameInstability))

        # Read the base flow
        print('Reading base flow: ')

        if os.path.isfile(os.path.join(caseFolder, 'baseflow.mat')):
            baseflow = scipy.io.loadmat(os.path.join(caseFolder, 'baseflow.mat'))
            print(f'{caseFolder}/baseflow.mat')
        else:
            nStep = checkPreviousRun(caseFolder)
            baseflow = scipy.io.loadmat(f'{caseFolder}/flow_{nStep:010d}.mat')
            print(f'{caseFolder}/flow_{nStep:010d}.mat')

        # If needed, replicate flow to make it 3D
        if mesh['z']['n'] > 1 and baseflow['U'].shape[2] == 1:
            baseflow['U'] = np.tile(baseflow['U'], (1, 1, mesh['z']['n']))
            baseflow['V'] = np.tile(baseflow['V'], (1, 1, mesh['z']['n']))
            baseflow['W'] = np.tile(baseflow['W'], (1, 1, mesh['z']['n']))
            baseflow['R'] = np.tile(baseflow['R'], (1, 1, mesh['z']['n']))
            baseflow['E'] = np.tile(baseflow['E'], (1, 1, mesh['z']['n']))

        # Get the physical flow region
        flowRegion = ~np.isnan(baseflow['U'])  # Get region outside of walls

        if removeBufferZone:
            flowRegion[np.r_[0:mesh['x']['buffer']['i']['n'], -mesh['x']['buffer']['f']['n']:], :, :] = False  # Remove buffer zones
            flowRegion[:, np.r_[0:mesh['y']['buffer']['i']['n'], -mesh['y']['buffer']['f']['n']:], :] = False
            flowRegion[:, :, np.r_[0:mesh['z']['buffer']['i']['n'], -mesh['z']['buffer']['f']['n']:]] = False

        if mesh['z']['n'] == 1:
            N = 4 * np.sum(flowRegion)
        elif not singleBeta:
            N = 5 * np.sum(flowRegion)
        else:
            N = 5 * np.sum(flowRegion[:, :, 0])

        # Allocate variables
        print('Allocating variables')
        zeta = np.zeros((N, M))
        H = np.zeros((M, M))

        # Check for previous results and reload data
        if os.path.isfile(os.path.join(caseFolder, f'results-{caseNameInstability}.mat')):
            print('Previous results file found, resuming')
            previousResults = scipy.io.loadmat(os.path.join(caseFolder, f'results-{caseNameInstability}.mat'))
            
            M0 = min(previousResults['M'][0][0], M)
            zeta[:, :M0] = previousResults['zeta'][:, :M0]
            H[:M0, :M0-1] = previousResults['H'][:M0, :M0-1]
        else:
            M0 = 1
            zeta[:, 0] = calcInitialDisturb(mesh['X'], mesh['Y'], mesh['Z'], singleBeta, flowRegion)

        # Set the initial SFD state
        if resumeSFD:
            t = 0
            U = baseflow['U']
            V = baseflow['V']
            W = baseflow['W']
            R = baseflow['R']
            E = baseflow['E']
            scipy.io.savemat(os.path.join(caseFolder, 'Instability', caseNameInstability, 'meanflowSFD.mat'), {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t})

        # Iterate
        if firstOrder:  # If using a first order approximation, compute the baseflow drift
            print('Iteration 0')
            t = 0
            
            U = baseflow['U']
            V = baseflow['V']
            W = baseflow['W']
            R = baseflow['R']
            E = baseflow['E']
            scipy.io.savemat(os.path.join(caseFolder, 'Instability', caseNameInstability, 'flow_0000000000.mat'), {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t})
            
            os.system(f'(cd {os.path.join(caseFolder, "Instability", caseNameInstability, "bin")} && mpirun -np {p_row * p_col} main {caseNameInstability}) {suppressOutput}')
            flowMinus = scipy.io.loadmat(os.path.join(caseFolder, 'Instability', caseNameInstability, 'flow_0000000001.mat'))
            os.remove(os.path.join(caseFolder, 'Instability', caseNameInstability, 'flow_0000000001.mat'))
            
            UM = flow2vec(flowMinus, flowRegion, singleBeta)

        elapsedTime = []
        for k in range(M0, M):
            print('\n')
            print(f'Iteration {k} of {M}')
            
            # Get current disturbance
            disturbance = vec2flow(zeta[:, k] * epsilon * np.sqrt(N), flowRegion, singleBeta)
            
            # Save flows, run DNS and read results
            t = 0
            
            U = baseflow['U'] + disturbance['U']
            V = baseflow['V'] + disturbance['V']
            W = baseflow['W'] + disturbance['W']
            R = baseflow['R'] + disturbance['R']
            E = baseflow['E'] + disturbance['E']
            scipy.io.savemat(os.path.join(caseFolder, 'Instability', caseNameInstability, 'flow_0000000000.mat'), {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t})
            
            os.system(f'(cd {os.path.join(caseFolder, "Instability", caseNameInstability, "bin")} && mpirun -np {p_row * p_col} main {caseNameInstability}) {suppressOutput}')
            flowPlus = scipy.io.loadmat(os.path.join(caseFolder, 'Instability', caseNameInstability, 'flow_0000000001.mat'))
            os.remove(os.path.join(caseFolder, 'Instability', caseNameInstability, 'flow_0000000001.mat'))
            
            UP = flow2vec(flowPlus, flowRegion, singleBeta)
            
            if not firstOrder:
                U = baseflow['U'] - disturbance['U']
                V = baseflow['V'] - disturbance['V']
                W = baseflow['W'] - disturbance['W']
                R = baseflow['R'] - disturbance['R']
                E = baseflow['E'] - disturbance['E']
                scipy.io.savemat(os.path.join(caseFolder, 'Instability', caseNameInstability, 'flow_0000000000.mat'), {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t})
                
                os.system(f'(cd {os.path.join(caseFolder, "Instability", caseNameInstability, "bin")} && mpirun -np {p_row * p_col} main {caseNameInstability}) {suppressOutput}')
                flowMinus = scipy.io.loadmat(os.path.join(caseFolder, 'Instability', caseNameInstability, 'flow_0000000001.mat'))
                os.remove(os.path.join(caseFolder, 'Instability', caseNameInstability, 'flow_0000000001.mat'))

                UM = flow2vec(flowMinus, flowRegion, singleBeta)

                # Compute matrices
                B = (UP - UM) / (2 * epsilon * np.sqrt(N))
            else:
                B = (UP - UM) / (epsilon * np.sqrt(N))
            
            uPrime = B

            for j in range(1, k + 1):
                H[j - 1, k] = np.dot(zeta[:, j - 1], B)
                uPrime -= H[j - 1, k] * zeta[:, j - 1]

            if k < M:
                H[k, k - 1] = np.linalg.norm(uPrime)
                zeta[:, k] = uPrime / H[k, k - 1]

            if np.any(np.isnan(H)) or np.any(np.isinf(H)):
                raise ValueError('DNS solution has failed')

            # Save results to file
            if k % saveEvery == 0 and k != M0:
                actualM = M
                actualH = H.copy()
                actualZeta = zeta.copy()
                
                M = k
                
                H[M + 1:, :] = []
                H[:, M + 1:] = []
                zeta[:, M + 1:] = []
                
                print(f'Saving results to results-{caseNameInstability}.mat')
                scipy.io.savemat(os.path.join(caseFolder, f'results-{caseNameInstability}.mat'), {'zeta': zeta, 'H': H, 'M': M})

                zeta = actualZeta
                H = actualH
                M = actualM

                del actualZeta, actualM

            # Compute modes
            if k in nKrilov:
                print('Computing modes')

                psiH, lambdaH = np.linalg.eig(H[:k, :k])
                lambdaB = np.diag(lambdaH)

                lambdaA = np.log(lambdaB) / dT
                index = np.argsort(np.real(lambdaA))[::-1]
                lambdaA = lambdaA[index]

                psiH = psiH[:, index]

                psiA = zeta[:, :k] @ psiH[:, :min(len(nModes), k)]
                
                modes = vec2modes(psiA, flowRegion, singleBeta, min(len(nModes), k), lambdaA, removeCC)
                modes['lambda'] = lambdaA
                modes['X'] = mesh['X']
                modes['Y'] = mesh['Y']
                modes['Z'] = mesh['Z']
                modes['nModes'] = nModes
                
                # Save modes
                print(f'Saving modes to modes-{caseNameInstability}-M{k}.mat')
                scipy.io.savemat(os.path.join(caseFolder, f'modes-{caseNameInstability}-M{k}.mat'), modes)

            # Print remaining time
            elapsedTime.append(time.time())
            remainingTime = (elapsedTime[-1] - elapsedTime[0]) / (len(elapsedTime) - 1) * (M - k)
            if len(elapsedTime) > 100:
                elapsedTime.pop(0)
            print(f'Elapsed time: {sec2str(elapsedTime[-1])}')
            print(f'Estimated time remaining: {sec2str(remainingTime)}')

            # Print leading eigenvalues
            if k % printEigsEvery == 0:
                leadingEigs = eigs(H, min(k, printEigsN), which='largestabs', tol=1e-4, maxiter=max(300, 2 * H.shape[0]), failure_treatment='keep')
                leadingEigs = np.log(leadingEigs) / dT
                print(f'Leading eigenvalues : {leadingEigs}')
                if 'leadingEigsPrev' in locals():
                    print(f'Change : {np.abs((leadingEigs[:len(leadingEigsPrev)] - leadingEigsPrev) / leadingEigs[:len(leadingEigsPrev)])}')
                leadingEigsPrev = leadingEigs


# Extra functions
def flow2vec(flow, flowRegion, singleBeta):
    threeD = flowRegion.shape[2] > 1
    
    if not threeD:  # 2D flow
        vec = np.concatenate((flow['U'][flowRegion], flow['V'][flowRegion], flow['R'][flowRegion], flow['E'][flowRegion]))
        
    elif not singleBeta:  # Fully 3D flow
        vec = np.concatenate((flow['U'][flowRegion], flow['V'][flowRegion], flow['W'][flowRegion], flow['R'][flowRegion], flow['E'][flowRegion]))
    
    else:  # Single beta, 3D flow
        if flow['U'].shape[2] > 1:
            nz = flowRegion.shape[2]
            Ztemp = np.linspace(0, 2 * np.pi, nz + 1)[:-1]
            Ztemp = np.transpose(Ztemp, (0, 2, 1))

            cosZ = np.cos(Ztemp) * 2 / nz
            sinZ = np.sin(Ztemp) * 2 / nz
            
            flow['U'] = np.sum(cosZ * flow['U'], axis=2)
            flow['V'] = np.sum(cosZ * flow['V'], axis=2)
            flow['W'] = np.sum(sinZ * flow['W'], axis=2)
            flow['R'] = np.sum(cosZ * flow['R'], axis=2)
            flow['E'] = np.sum(cosZ * flow['E'], axis=2)
        
        flowRegion = flowRegion[:, :, 0]
        
        vec = np.concatenate((flow['U'][flowRegion], flow['V'][flowRegion], flow['W'][flowRegion], flow['R'][flowRegion], flow['E'][flowRegion]))
    
    return vec

def vec2flow(vec, flowRegion, singleBeta):
    nx, ny, nz = flowRegion.shape
    N = np.sum(flowRegion)
    
    threeD = nz > 1
    
    if not threeD:  # 2D flow
        flow = {
            'U': np.zeros((nx, ny)),
            'V': np.zeros((nx, ny)),
            'W': np.zeros((nx, ny)),
            'R': np.zeros((nx, ny)),
            'E': np.zeros((nx, ny))
        }
        
        flow['U'][flowRegion] = vec[:N]
        flow['V'][flowRegion] = vec[N:2*N]
        flow['R'][flowRegion] = vec[2*N:3*N]
        flow['E'][flowRegion] = vec[3*N:4*N]
        
    elif not singleBeta:  # Fully 3D flow
        flow = {
            'U': np.zeros((nx, ny, nz)),
            'V': np.zeros((nx, ny, nz)),
            'W': np.zeros((nx, ny, nz)),
            'R': np.zeros((nx, ny, nz)),
            'E': np.zeros((nx, ny, nz))
        }
        
        flow['U'][flowRegion] = vec[:N]
        flow['V'][flowRegion] = vec[N:2*N]
        flow['W'][flowRegion] = vec[2*N:3*N]
        flow['R'][flowRegion] = vec[3*N:4*N]
        flow['E'][flowRegion] = vec[4*N:5*N]
    
    else:  # Single beta, 3D flow
        flow = {
            'U': np.zeros((nx, ny)),
            'V': np.zeros((nx, ny)),
            'W': np.zeros((nx, ny)),
            'R': np.zeros((nx, ny)),
            'E': np.zeros((nx, ny))
        }
        
        N = round(N / nz)
        
        flow['U'][flowRegion[:, :, 0]] = vec[:N]
        flow['V'][flowRegion[:, :, 0]] = vec[N:2*N]
        flow['W'][flowRegion[:, :, 0]] = vec[2*N:3*N]
        flow['R'][flowRegion[:, :, 0]] = vec[3*N:4*N]
        flow['E'][flowRegion[:, :, 0]] = vec[4*N:5*N]
    
        Ztemp = np.linspace(0, 2 * np.pi, nz + 1)[:-1]
        Ztemp = np.transpose(Ztemp, (0, 2, 1))
        
        cosZ = np.cos(Ztemp)
        sinZ = np.sin(Ztemp)
        
        flow['U'] *= cosZ
        flow['V'] *= cosZ
        flow['W'] *= sinZ
        flow['R'] *= cosZ
        flow['E'] *= cosZ
    
    return flow

def vec2modes(vec, flowRegion, singleBeta, nModes, lambda_, removeCC):
    nx, ny, nz = flowRegion.shape
    N = np.sum(flowRegion)
    
    threeD = nz > 1
    
    if not threeD:  # 2D flow
        modes = {
            'U': np.full((nx, ny, 1, len(nModes)), np.nan),
            'V': np.full((nx, ny, 1, len(nModes)), np.nan),
            'W': np.full((nx, ny, 1, len(nModes)), np.nan),
            'R': np.full((nx, ny, 1, len(nModes)), np.nan),
            'E': np.full((nx, ny, 1, len(nModes)), np.nan)
        }
        
        U = np.full((nx, ny), np.nan)
        V = np.full((nx, ny), np.nan)
        R = np.full((nx, ny), np.nan)
        E = np.full((nx, ny), np.nan)
        
        for i in range(len(nModes)):
            if np.imag(lambda_[i]) >= 0:
                vec[:, nModes[i]] /= np.max(np.abs(vec[:, nModes[i]]))
                
                U[flowRegion] = vec[:N, nModes[i]]
                V[flowRegion] = vec[N:2*N, nModes[i]]
                R[flowRegion] = vec[2*N:3*N, nModes[i]]
                E[flowRegion] = vec[3*N:4*N, nModes[i]]
                
                modes['U'][:, :, 0, i] = U
                modes['V'][:, :, 0, i] = V
                modes['R'][:, :, 0, i] = R
                modes['E'][:, :, 0, i] = E
        
    elif not singleBeta:  # Fully 3D flow
        modes = {
            'U': np.full((nx, ny, nz, len(nModes)), np.nan),
            'V': np.full((nx, ny, nz, len(nModes)), np.nan),
            'W': np.full((nx, ny, nz, len(nModes)), np.nan),
            'R': np.full((nx, ny, nz, len(nModes)), np.nan),
            'E': np.full((nx, ny, nz, len(nModes)), np.nan)
        }
        
        U = np.full((nx, ny, nz), np.nan)
        V = np.full((nx, ny, nz), np.nan)
        W = np.full((nx, ny, nz), np.nan)
        R = np.full((nx, ny, nz), np.nan)
        E = np.full((nx, ny, nz), np.nan)
        
        for i in range(len(nModes)):
            if np.imag(lambda_[i]) >= 0:
                vec[:, nModes[i]] /= np.max(np.abs(vec[:, nModes[i]]))
                
                U[flowRegion] = vec[:N, nModes[i]]
                V[flowRegion] = vec[N:2*N, nModes[i]]
                W[flowRegion] = vec[2*N:3*N, nModes[i]]
                R[flowRegion] = vec[3*N:4*N, nModes[i]]
                E[flowRegion] = vec[4*N:5*N, nModes[i]]
                
                modes['U'][:, :, :, i] = U
                modes['V'][:, :, :, i] = V
                modes['W'][:, :, :, i] = W
                modes['R'][:, :, :, i] = R
                modes['E'][:, :, :, i] = E
        
    else:  # Single beta, 3D flow
        N = round(N / nz)
        flowRegion = flowRegion[:, :, 0]
        
        modes = {
            'U': np.full((nx, ny, 1, len(nModes)), np.nan),
            'V': np.full((nx, ny, 1, len(nModes)), np.nan),
            'W': np.full((nx, ny, 1, len(nModes)), np.nan),
            'R': np.full((nx, ny, 1, len(nModes)), np.nan),
            'E': np.full((nx, ny, 1, len(nModes)), np.nan)
        }
        
        U = np.full((nx, ny), np.nan)
        V = np.full((nx, ny), np.nan)
        W = np.full((nx, ny), np.nan)
        R = np.full((nx, ny), np.nan)
        E = np.full((nx, ny), np.nan)
        
        for i in range(len(nModes)):
            if np.imag(lambda_[i]) >= 0:
                vec[:, nModes[i]] /= np.max(np.abs(vec[:, nModes[i]]))
                
                U[flowRegion] = vec[:N, nModes[i]]
                V[flowRegion] = vec[N:2*N, nModes[i]]
                W[flowRegion] = vec[2*N:3*N, nModes[i]]
                R[flowRegion] = vec[3*N:4*N, nModes[i]]
                E[flowRegion] = vec[4*N:5*N, nModes[i]]
                
                modes['U'][:, :, 0, i] = U
                modes['V'][:, :, 0, i] = V
                modes['W'][:, :, 0, i] = W
                modes['R'][:, :, 0, i] = R
                modes['E'][:, :, 0, i] = E
    
    return modes

def calcInitialDisturb(X, Y, Z, singleBeta, flowRegion):
    x0i = np.argmin(np.gradient(X))
    y0i = np.argmin(np.gradient(Y))
    z0i = np.argmin(np.gradient(Z))

    x0 = X[x0i]
    y0 = Y[y0i]
    z0 = Z[z0i]
    
    alphaX = 10 / (X[-1] - X[0]) ** 2
    alphaY = 10 / (Y[-1] - Y[0]) ** 2
    alphaZ = 10 / (Z[-1] - Z[0]) ** 2
    
    if len(Z) == 1:
        eta = np.exp(-(alphaX * (X[:, None] - x0) ** 2 + alphaY * (Y - y0) ** 2))
        
        flow = {
            'U': eta,
            'V': eta,
            'W': eta,
            'R': eta,
            'E': eta
        }
        
        zeta = flow2vec(flow, flowRegion, singleBeta)
        
    elif not singleBeta:
        eta = np.exp(-(alphaX * (X[:, None] - x0) ** 2 + alphaY * (Y - y0) ** 2 + alphaZ * (np.transpose(Z) - z0) ** 2))
        
        flow = {
            'U': eta,
            'V': eta,
            'W': eta,
            'R': eta,
            'E': eta
        }
        
        zeta = flow2vec(flow, flowRegion, singleBeta)
        
    else:
        eta = np.exp(-(alphaX * (X[:, None] - x0) ** 2 + alphaY * (Y - y0) ** 2))
        
        flow = {
            'U': eta,
            'V': eta,
            'W': eta,
            'R': eta,
            'E': eta
        }
        
        zeta = flow2vec(flow, flowRegion, singleBeta)
    
    zeta = zeta / np.linalg.norm(zeta)
    
    return zeta

def sec2str(sec):
    sec = round(sec)
    if sec < 60:
        return f'{sec}s'
    else:
        min_ = sec // 60
        sec = sec % 60
    if min_ < 60:
        return f'{min_}min {sec:02d}s'
    else:
        hour = min_ // 60
        min_ = min_ % 60
    if hour < 24:
        return f'{hour}h {min_:02d}min {sec:02d}s'
    else:
        day = hour // 24
        hour = hour % 24
    return f'{day}d {hour:02d}h {min_:02d}min {sec:02d}s'

def unpackStruct(structure):
    varList = []
    varNames = list(structure.keys())
    for varName in varNames:
        if isinstance(structure[varName], dict):
            varList2 = unpackStruct(structure[varName])
            for item in varList2:
                varList.append(f'{varName}.{item}')
        else:
            varList.append(varName)
    return varList

def cleanUpFunction(caseFolder, caseNameInstability=None):
    import os
    if caseNameInstability is None:  # Restore the original files
        os.system(f'mv {caseFolder}/mesh.mat.backup {caseFolder}/mesh.mat >/dev/null 2>&1')
        os.system(f'rm -r {caseFolder}/bin >/dev/null 2>&1')
        os.system(f'mv {caseFolder}/bin_backup {caseFolder}/bin >/dev/null 2>&1')
    else:  # Clean up instability folder
        os.system(f'rm -r {caseFolder}/Instability/{caseNameInstability}')
        temp = os.listdir(f'{caseFolder}/Instability')
        if len(temp) == 2:  # Check if Instability dir is empty and erase it
            os.system(f'rm -r "$(readlink -f {caseFolder}/Instability)"')
            os.system(f'rm -r {caseFolder}/Instability')

def get_p_row(baseFile, p_row_max):
    import importlib
    
    importlib.import_module('source')
    importlib.import_module('source.boundaries')
    
    eval(baseFile)
    p_row = p_row_max
    p_col = 1
    while p_row > 0:
        try:
            p_row
            meshAddFixedPoints()

            # Run mesh generator
            mesh = {}
            mesh['X'], mesh['x'], mesh['nx'] = generateMesh(domain['xi'], domain['xf'], mesh['x'], 'X')
            mesh['Y'], mesh['y'], mesh['ny'] = generateMesh(domain['yi'], domain['yf'], mesh['y'], 'Y')
            mesh['Z'], mesh['z'], mesh['nz'] = generateMesh(domain['zi'], domain['zf'], mesh['z'], 'Z')

            # Select boundary conditions
            boundary, mesh = getBoundaryConditions(flowType, mesh, flowParameters, [numMethods['neumannOrder'], numMethods['neumann2Order']])

            domainSlicesY = getDomainSlices(mesh['ny'], p_row)
            domainSlicesZ = getDomainSlices(mesh['nz'], p_col)

            initBoundaries(boundary, mesh, domainSlicesY, domainSlicesZ, p_row, p_col)
            
            return
        except:
            p_row -= 1



