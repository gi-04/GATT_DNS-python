import numpy as np
from scipy.interpolate import interp1d, griddata
from scipy.integrate import solve_ivp

def generate_initial_flow(mesh, flow_parameters, initial_flow, walls, flow_name):
    # This function defines each type of initial flow. New types may be added here if needed

    nx = mesh['nx']
    ny = mesh['ny']
    nz = mesh['nz']

    X = mesh['X']
    Y = mesh['Y']
    Z = mesh['Z']

    gamma = flow_parameters['gamma']
    Ma = flow_parameters['Ma']
    Re = flow_parameters['Re']

    E0 = 1 / ((gamma**2 - gamma) * Ma**2)

    if initial_flow['type'] == 'uniform':
        
        U = np.ones((nx, ny, nz))
        V = np.zeros((nx, ny, nz))
        W = np.zeros((nx, ny, nz))
        R = np.ones((nx, ny, nz))
        E = np.ones((nx, ny, nz)) * E0
        
        y0ind = np.argmin(np.abs(Y))
        
        if len(initial_flow) == 4 or 'boundaryLayer' in flow_name:
            U[:, :y0ind, :] = 0
        
        if 'U0' in initial_flow:
            U = initial_flow['U0'] * U
        
    elif initial_flow['type'] == 'poiseuille':
        
        U = np.zeros((nx, ny, nz))
        V = np.zeros((nx, ny, nz))
        W = np.zeros((nx, ny, nz))
        R = np.ones((nx, ny, nz))
        E = np.ones((nx, ny, nz)) * E0

        eta = (Y - Y[0]) / (Y[-1] - Y[0])

        u0 = flow_parameters['U0']
        u1 = flow_parameters['lowerWallVelocity']
        u2 = flow_parameters['upperWallVelocity']

        U += (-6 * u0 + 3 * u1 + 3 * u2) * eta**2 + (6 * u0 - 4 * u1 - 2 * u2) * eta + u1

    elif initial_flow['type'] == 'blasius':

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
        if 'blasiusFit' in initial_flow:
            walls2D = np.any(walls, axis=2)
            for i in range(nx):
                Y0[i] = Y[np.where(walls2D[i, :] == 0)[0][0]]
            for i in range(1, nx):
                dx = X[i] - X[i - 1]
                if Y0[i - 1] - Y0[i] > dx * initial_flow['blasiusFit']:
                    Y0[i] = Y0[i - 1] - dx * initial_flow['blasiusFit']
            for i in range(nx - 1, 0, -1):
                dx = X[i + 1] - X[i]
                if Y0[i + 1] - Y0[i] > dx * initial_flow['blasiusFit']:
                    Y0[i] = Y0[i + 1] - dx * initial_flow['blasiusFit']

        for i in range(nx):
            xi = X[i]
            if xi > 0:
                theta = 0.664 * np.sqrt(xi / Re)
                U[i, :] = interp1d(ybl * theta / thetabl + Y0[i], ubl, kind='cubic', fill_value="extrapolate")(Y)
                U[i, Y < Y0[i]] = 0

        U[U > 1] = 1
        U[np.isnan(U)] = 1

        U = np.repeat(U[:, :, np.newaxis], nz, axis=2)
        
    elif initial_flow['type'] == 'compressibleBL_isothermal':
        compressible_bl_flow = calc_compressible_bl(flow_parameters, False, mesh)

        U = compressible_bl_flow['U']
        V = compressible_bl_flow['V']
        W = compressible_bl_flow['W']
        E = compressible_bl_flow['E']
        R = compressible_bl_flow['R']

    elif initial_flow['type'] == 'compressibleBL_adiabatic':
        compressible_bl_flow = calc_compressible_bl(flow_parameters, True, mesh)

        U = compressible_bl_flow['U']
        V = compressible_bl_flow['V']
        W = compressible_bl_flow['W']
        E = compressible_bl_flow['E']
        R = compressible_bl_flow['R']

    elif initial_flow['type'] == 'file':
        
        # If a folder is given, find the last flow file inside it
        if initial_flow['flowFile'][-1] == '/':
            nStep = check_previous_run(initial_flow['flowFile'][:-1])
            if nStep:
                initial_flow['flowFile'] = f"{initial_flow['flowFile']}flow_{nStep:010d}.mat"
            else:
                initial_flow['flowFile'] = f"{initial_flow['flowFile']}baseflow.mat"

        flow_file = loadmat(initial_flow['flowFile'])
        
        # If a mesh file is given, interpolate from it
        if 'meshFile' in initial_flow:

            # If a folder is given, use the mesh file inside it
            if initial_flow['meshFile'][-1] == '/':
                initial_flow['meshFile'] += 'mesh.mat'

            mesh_file = loadmat(initial_flow['meshFile'])
            Xfile = mesh_file['X']
            Yfile = mesh_file['Y']
            Zfile = mesh_file['Z']
            
            Ufile = flow_file['U']
            Vfile = flow_file['V']
            Wfile = flow_file['W']
            Rfile = flow_file['R']
            Efile = flow_file['E']
            
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
            if 'changeMach' in initial_flow and initial_flow['changeMach']:
                if mesh_file['flowParameters']['Ma'] != Ma:
                    print(f'Initial flow file has a Mach number of {mesh_file["flowParameters"]["Ma"]} which will be changed to {Ma} to match the current simulation')
                    E *= mesh_file['flowParameters']['Ma']**2 / Ma**2

        else:
            # If no mesh file is given, just load the values and check the size
            U = flow_file['U']
            if nx != U.shape[0] or ny != U.shape[1] or (nz != U.shape[2] and U.shape[2] != 1):
                raise ValueError(f'Mesh size in initial flow file is not consistent with current mesh ({U.shape[0]}x{U.shape[1]}x{U.shape[2]} -> {nx}x{ny}x{nz})\nConsider changing the parameters file or providing a mesh file')
            V = flow_file['V']
            W = flow_file['W']
            R = flow_file['R']
            E = flow_file['E']
            
            # If the mesh provided was 2D and the case is 3D, replicate it
            if nz == 1 and U.shape[2]:
                U = np.repeat(U[:, :, np.newaxis], nz, axis=2)
                V = np.repeat(V[:, :, np.newaxis], nz, axis=2)
                W = np.repeat(W[:, :, np.newaxis], nz, axis=2)
                R = np.repeat(R[:, :, np.newaxis], nz, axis=2)
                E = np.repeat(E[:, :, np.newaxis], nz, axis=2)
        
    if 'addNoise' in initial_flow:

        if 'noiseType' not in initial_flow or initial_flow['noiseType'] == 'rand':
            noiseU = initial_flow['addNoise'] * np.random.randn(nx, ny, nz)
            noiseV = initial_flow['addNoise'] * np.random.randn(nx, ny, nz)
            noiseW = initial_flow['addNoise'] * np.random.randn(nx, ny, nz)
            noiseR = initial_flow['addNoise'] * np.random.randn(nx, ny, nz)
            noiseE = initial_flow['addNoise'] * np.random.randn(nx, ny, nz)
        elif initial_flow['noiseType'] == 'uniform':
            noiseU = initial_flow['addNoise'] * np.ones((nx, ny, nz))
            noiseV = initial_flow['addNoise'] * np.ones((nx, ny, nz))
            noiseW = initial_flow['addNoise'] * np.ones((nx, ny, nz))
            noiseR = initial_flow['addNoise'] * np.ones((nx, ny, nz))
            noiseE = initial_flow['addNoise'] * np.ones((nx, ny, nz))

        if 'noiseCenter' in initial_flow:
            if 'noiseSigma' not in initial_flow:
                initial_flow['noiseSigma'] = [(X[-1] - X[0])**2 / 10, (Y[-1] - Y[0])**2 / 10, (Z[-1] - Z[0])**2 / 10]
            x0 = initial_flow['noiseCenter'][0]
            y0 = initial_flow['noiseCenter'][1]
            z0 = initial_flow['noiseCenter'][2]
            
            sigmaX = initial_flow['noiseSigma'][0]
            sigmaY = initial_flow['noiseSigma'][1]
            sigmaZ = initial_flow['noiseSigma'][2]
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



