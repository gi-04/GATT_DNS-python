# Required python packages
import numpy as np
from scipy.interpolate import interp1d, griddata
from scipy.integrate import solve_ivp
from scipy.io import loadmat
# ADditional functions and files
from calcCompressibleBL import calcCompressibleBL
from checkPreviousRun import checkPreviousRun

def generateInitialFlow(mesh, flowParameters, initialFlow, walls, flow_name):
    # This function defines each type of initial flow. New types may be added here if needed

    nx = mesh['nx']
    ny = mesh['ny']
    nz = mesh['nz']

    X = mesh['X']
    Y = mesh['Y']
    Z = mesh['Z']

    gamma = flowParameters['gamma']
    Ma = flowParameters['Ma']
    Re = flowParameters['Re']

    E0 = 1 / ((gamma**2 - gamma) * Ma**2)

    if initialFlow['type'] == 'uniform':
        
        U = np.ones((nx, ny, nz))
        V = np.zeros((nx, ny, nz))
        W = np.zeros((nx, ny, nz))
        R = np.ones((nx, ny, nz))
        E = np.ones((nx, ny, nz)) * E0
        
        y0ind = np.argmin(np.abs(Y))
        
        if len(initialFlow) == 4 or 'boundaryLayer' in flow_name:
            U[:, :y0ind, :] = 0
        
        if 'U0' in initialFlow:
            U = initialFlow['U0'] * U
        
    elif initialFlow['type'] == 'poiseuille':
        
        U = np.zeros((nx, ny, nz))
        V = np.zeros((nx, ny, nz))
        W = np.zeros((nx, ny, nz))
        R = np.ones((nx, ny, nz))
        E = np.ones((nx, ny, nz)) * E0

        eta = (Y - Y[0]) / (Y[-1] - Y[0])

        u0 = flowParameters['U0']
        u1 = flowParameters['lowerWallVelocity']
        u2 = flowParameters['upperWallVelocity']

        U += (-6 * u0 + 3 * u1 + 3 * u2) * eta**2 + (6 * u0 - 4 * u1 - 2 * u2) * eta + u1

    elif initialFlow['type'] == 'blasius':

        ybl = np.arange(0, 10.0001, 0.0001)
        ubl = blasius(ybl)
        thetabl = 0.664155332943009

        e0 = 1 / ((gamma**2 - gamma) * Ma**2)

        nx = len(X)
        ny = len(Y)
        nz = len(Z)

        R = np.ones((nx, ny, nz))
        U = np.ones((nx, ny))
        V = np.zeros((nx, ny, nz))
        W = np.zeros((nx, ny, nz))
        E = np.ones((nx, ny, nz)) * e0

        # If blasiusFit defined, the blasius boundary layer will fit to the wall.
        # The value of blasiusFit defines the maximum slope of Y0 for the blasius profile

        Y0 = np.zeros(nx)
        if 'blasiusFit' in initialFlow:
            walls2D = np.any(walls, axis=2)
            for i in range(nx):
                Y0[i] = Y[np.where(walls2D[i, :] == 0)[0][0]]
            for i in range(1, nx):
                dx = X[i] - X[i - 1]
                if Y0[i - 1] - Y0[i] > dx * initialFlow['blasiusFit']:
                    Y0[i] = Y0[i - 1] - dx * initialFlow['blasiusFit']
            for i in range(nx - 1, 0, -1):
                dx = X[i + 1] - X[i]
                if Y0[i + 1] - Y0[i] > dx * initialFlow['blasiusFit']:
                    Y0[i] = Y0[i + 1] - dx * initialFlow['blasiusFit']

        for i in range(nx):
            xi = X[i]
            if xi > 0:
                theta = 0.664 * np.sqrt(xi / Re)
                U[i, :] = interp1d(ybl * theta / thetabl + Y0[i], ubl, kind='cubic', fill_value="extrapolate")(Y)
                U[i, Y < Y0[i]] = 0

        U[U > 1] = 1
        U[np.isnan(U)] = 1

        U = np.repeat(U[:, :, np.newaxis], nz, axis=2)
        
    elif initialFlow['type'] == 'compressibleBL_isothermal':
        compressibleBL_flow = calcCompressibleBL(flowParameters, False, mesh)

        U = compressibleBL_flow['U']
        V = compressibleBL_flow['V']
        W = compressibleBL_flow['W']
        E = compressibleBL_flow['E']
        R = compressibleBL_flow['R']

    elif initialFlow['type'] == 'compressibleBL_adiabatic':
        compressibleBL_flow = calcCompressibleBL(flowParameters, True, mesh)

        U = compressibleBL_flow['U']
        V = compressibleBL_flow['V']
        W = compressibleBL_flow['W']
        E = compressibleBL_flow['E']
        R = compressibleBL_flow['R']

    elif initialFlow['type'] == 'file':
        
        # If a folder is given, find the last flow file inside it
        if initialFlow['flowFile'][-1] == '/':
            nStep = checkPreviousRun(initialFlow['flowFile'][:-1])
            if nStep:
                initialFlow['flowFile'] = f"{initialFlow['flowFile']}flow_{nStep:010d}.mat"
            else:
                initialFlow['flowFile'] = f"{initialFlow['flowFile']}baseflow.mat"

        flowFile = loadmat(initialFlow['flowFile'])
        
        # If a mesh file is given, interpolate from it
        if 'meshFile' in initialFlow:

            # If a folder is given, use the mesh file inside it
            if initialFlow['meshFile'][-1] == '/':
                initialFlow['meshFile'] += 'mesh.mat'

            meshFile = loadmat(initialFlow['meshFile'])
            Xfile = meshFile['X']
            Yfile = meshFile['Y']
            Zfile = meshFile['Z']
            
            Ufile = flowFile['U']
            Vfile = flowFile['V']
            Wfile = flowFile['W']
            Rfile = flowFile['R']
            Efile = flowFile['E']
            
            # Change NaNs to default values
            Ufile[np.isnan(Ufile)] = 0
            Vfile[np.isnan(Vfile)] = 0
            Wfile[np.isnan(Wfile)] = 0
            Rfile[np.isnan(Rfile)] = 1
            Efile[np.isnan(Efile)] = E0
            
            # Create a temporary mesh and crop it to avoid extrapolations
            Xmesh, Ymesh, Zmesh = np.meshgrid(X, Y, Z, indexing='ij')
            Xmesh[Xmesh < Xfile[0]] = Xfile[0]
            Xmesh[Xmesh > Xfile[-1]] = Xfile[-1]
            Ymesh[Ymesh < Yfile[0]] = Yfile[0]
            Ymesh[Ymesh > Yfile[-1]] = Yfile[-1]
            Zmesh[Zmesh < Zfile[0]] = Zfile[0]
            Zmesh[Zmesh > Zfile[-1]] = Zfile[-1]
            
            # Interpolate
            if nz == 1 or len(Zfile) == 1:  # In a 2D case
                U = griddata((Xfile.flatten(), Yfile.flatten()), Ufile.flatten(), (Xmesh, Ymesh), method='linear')
                V = griddata((Xfile.flatten(), Yfile.flatten()), Vfile.flatten(), (Xmesh, Ymesh), method='linear')
                W = griddata((Xfile.flatten(), Yfile.flatten()), Wfile.flatten(), (Xmesh, Ymesh), method='linear')
                R = griddata((Xfile.flatten(), Yfile.flatten()), Rfile.flatten(), (Xmesh, Ymesh), method='linear')
                E = griddata((Xfile.flatten(), Yfile.flatten()), Efile.flatten(), (Xmesh, Ymesh), method='linear')
            else:  # In a 3D case
                U = griddata((Xfile.flatten(), Yfile.flatten(), Zfile.flatten()), Ufile.flatten(), (Xmesh, Ymesh, Zmesh), method='linear')
                V = griddata((Xfile.flatten(), Yfile.flatten(), Zfile.flatten()), Vfile.flatten(), (Xmesh, Ymesh, Zmesh), method='linear')
                W = griddata((Xfile.flatten(), Yfile.flatten(), Zfile.flatten()), Wfile.flatten(), (Xmesh, Ymesh, Zmesh), method='linear')
                R = griddata((Xfile.flatten(), Yfile.flatten(), Zfile.flatten()), Rfile.flatten(), (Xmesh, Ymesh, Zmesh), method='linear')
                E = griddata((Xfile.flatten(), Yfile.flatten(), Zfile.flatten()), Efile.flatten(), (Xmesh, Ymesh, Zmesh), method='linear')
            
            # Change Mach number if needed
            if 'changeMach' in initialFlow and initialFlow['changeMach']:
                if meshFile['flowParameters']['Ma'] != Ma:
                    print(f'Initial flow file has a Mach number of {meshFile["flowParameters"]["Ma"]} which will be changed to {Ma} to match the current simulation')
                    E *= meshFile['flowParameters']['Ma']**2 / Ma**2

        else:
            # If no mesh file is given, just load the values and check the size
            U = flowFile['U']
            if nx != U.shape[0] or ny != U.shape[1] or (nz != U.shape[2] and U.shape[2] != 1):
                raise ValueError(f'Mesh size in initial flow file is not consistent with current mesh ({U.shape[0]}x{U.shape[1]}x{U.shape[2]} -> {nx}x{ny}x{nz})\nConsider changing the parameters file or providing a mesh file')
            V = flowFile['V']
            W = flowFile['W']
            R = flowFile['R']
            E = flowFile['E']
            
            # If the mesh provided was 2D and the case is 3D, replicate it
            if nz == 1 and U.shape[2]:
                U = np.repeat(U[:, :, np.newaxis], nz, axis=2)
                V = np.repeat(V[:, :, np.newaxis], nz, axis=2)
                W = np.repeat(W[:, :, np.newaxis], nz, axis=2)
                R = np.repeat(R[:, :, np.newaxis], nz, axis=2)
                E = np.repeat(E[:, :, np.newaxis], nz, axis=2)
        
    if 'addNoise' in initialFlow:

        if 'noiseType' not in initialFlow or initialFlow['noiseType'] == 'rand':
            noiseU = initialFlow['addNoise'] * np.random.randn(nx, ny, nz)
            noiseV = initialFlow['addNoise'] * np.random.randn(nx, ny, nz)
            noiseW = initialFlow['addNoise'] * np.random.randn(nx, ny, nz)
            noiseR = initialFlow['addNoise'] * np.random.randn(nx, ny, nz)
            noiseE = initialFlow['addNoise'] * np.random.randn(nx, ny, nz)
        elif initialFlow['noiseType'] == 'uniform':
            noiseU = initialFlow['addNoise'] * np.ones((nx, ny, nz))
            noiseV = initialFlow['addNoise'] * np.ones((nx, ny, nz))
            noiseW = initialFlow['addNoise'] * np.ones((nx, ny, nz))
            noiseR = initialFlow['addNoise'] * np.ones((nx, ny, nz))
            noiseE = initialFlow['addNoise'] * np.ones((nx, ny, nz))

        if 'noiseCenter' in initialFlow:
            if 'noiseSigma' not in initialFlow:
                initialFlow['noiseSigma'] = [(X[-1] - X[0])**2 / 10, (Y[-1] - Y[0])**2 / 10, (Z[-1] - Z[0])**2 / 10]
            x0 = initialFlow['noiseCenter'][0]
            y0 = initialFlow['noiseCenter'][1]
            z0 = initialFlow['noiseCenter'][2]
            
            sigmaX = initialFlow['noiseSigma'][0]
            sigmaY = initialFlow['noiseSigma'][1]
            sigmaZ = initialFlow['noiseSigma'][2]
        else:
            x0 = 0
            y0 = 0
            z0 = 0
            sigmaX = np.inf
            sigmaY = np.inf
            sigmaZ = np.inf
        
        if nz == 1:
            sigmaZ = np.inf
        radius = (np.add.outer((X - x0)**2 / sigmaX, (Y - y0)**2 / sigmaY) +
                  np.add.outer((Z - z0)**2 / sigmaZ, np.zeros((nx, ny))))
        noiseGaussian = np.exp(-radius)
        
        U += noiseU * noiseGaussian
        V += noiseV * noiseGaussian
        if nz > 1:
            W += noiseW * noiseGaussian
        R += noiseR * noiseGaussian
        E += noiseE * noiseGaussian

    flow = {'U': U, 'V': V, 'W': W, 'R': R, 'E': E}

    return flow

def blasius(y):
    f0 = [0, 0, 0.33204312]
    sol = solve_ivp(blasius_eq, [y[0], y[-1]], f0, t_eval=y)
    return sol.y[1, :]

def blasius_eq(t, f):
    d = np.zeros(3)
    d[2] = -0.5 * f[0] * f[2]
    d[1] = f[2]
    d[0] = f[1]
    return d



