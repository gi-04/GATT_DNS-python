import numpy as np
from scipy.interpolate import interpn, interp1d
from scipy.integrate import odeint

def generate_initial_flow(mesh, flow_parameters, initial_flow, walls):
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

        U[:, :y0ind, :] = 0

        if 'U0' in initial_flow:
            U = initial_flow['U0'] * U

    elif initial_flow['type'] == 'poiseulle':
        U = np.zeros((nx, ny, nz))
        V = np.zeros((nx, ny, nz))
        W = np.zeros((nx, ny, nz))
        R = np.ones((nx, ny, nz))
        E = np.ones((nx, ny, nz)) * E0

        U += mesh['Y'] / (mesh['Y'][-1] - mesh['Y'][0])

        if 'U0' in initial_flow:
            U = initial_flow['U0'] * U

    elif initial_flow['type'] == 'blasius':
        ybl = np.arange(0, 10, 0.0001)
        ubl = blasius(ybl)
        thetabl = 0.664155332943009

        e0 = 1 / ((gamma**2 - gamma) * Ma**2)

        R = np.ones((nx, ny, nz))
        U = np.ones((nx, ny))
        V = np.zeros((nx, ny, nz))
        W = np.zeros((nx, ny, nz))
        E = np.ones((nx, ny, nz)) * e0

        Y0 = np.zeros(nx)
        if 'blasiusFit' in initial_flow:
            walls2D = np.any(walls, axis=2)
            for i in range(nx):
                Y0[i] = Y[np.where(walls2D[i, :] == 0)[0][0]]
            for i in range(1, nx):
                dx = X[i] - X[i - 1]
                if Y0[i - 1] - Y0[i] > dx * initial_flow['blasiusFit']:
                    Y0[i] = Y0[i - 1] - dx * initial_flow['blasiusFit']
            for i in range(nx - 2, -1, -1):
                dx = X[i + 1] - X[i]
                if Y0[i + 1] - Y0[i] > dx * initial_flow['blasiusFit']:
                    Y0[i] = Y0[i + 1] - dx * initial_flow['blasiusFit']

        for i in range(nx):
            xi = X[i]
            if xi > 0:
                theta = 0.664 * np.sqrt(xi / Re)
                U[i, :] = interp1d(ybl * theta / thetabl + Y0[i], ubl, kind='spline', fill_value="extrapolate")(Y)
                U[i, Y < Y0[i]] = 0

        U[U > 1] = 1
        U[np.isnan(U)] = 1

        U = np.tile(U[:, :, np.newaxis], (1, 1, nz))

    elif initial_flow['type'] == 'file':
        if initial_flow['flowFile'][-1] == '/':
            nStep = check_previous_run(initial_flow['flowFile'][:-1])
            if nStep is not None:
                initial_flow['flowFile'] = f"{initial_flow['flowFile']}flow_{nStep:010d}.mat"
            else:
                initial_flow['flowFile'] = f"{initial_flow['flowFile']}baseflow.mat"

        flow_file = np.load(initial_flow['flowFile'], allow_pickle=True).item()

        if 'meshFile' in initial_flow:
            if initial_flow['meshFile'][-1] == '/':
                initial_flow['meshFile'] = f"{initial_flow['meshFile']}mesh.mat"

            mesh_file = np.load(initial_flow['meshFile'], allow_pickle=True).item()
            Xfile = mesh_file['X']
            Yfile = mesh_file['Y']
            Zfile = mesh_file['Z']

            Ufile = flow_file['U']
            Vfile = flow_file['V']
            Wfile = flow_file['W']
            Rfile = flow_file['R']
            Efile = flow_file['E']

            Ufile[np.isnan(Ufile)] = 0
            Vfile[np.isnan(Vfile)] = 0
            Wfile[np.isnan(Wfile)] = 0
            Rfile[np.isnan(Rfile)] = 1
            Efile[np.isnan(Efile)] = E0

            Xmesh, Ymesh, Zmesh = np.meshgrid(X, Y, Z, indexing='ij')
            Xmesh[Xmesh < Xfile[0]] = Xfile[0]
            Xmesh[Xmesh > Xfile[-1]] = Xfile[-1]
            Ymesh[Ymesh < Yfile[0]] = Yfile[0]
            Ymesh[Ymesh > Yfile[-1]] = Yfile[-1]
            Zmesh[Zmesh < Zfile[0]] = Zfile[0]
            Zmesh[Zmesh > Zfile[-1]] = Zfile[-1]

            if nz == 1 or len(Zfile) == 1:
                U = interpn((Xfile, Yfile), Ufile, (Xmesh, Ymesh), method='linear')
                V = interpn((Xfile, Yfile), Vfile, (Xmesh, Ymesh), method='linear')
                W = interpn((Xfile, Yfile), Wfile, (Xmesh, Ymesh), method='linear')
                R = interpn((Xfile, Yfile), Rfile, (Xmesh, Ymesh), method='linear')
                E = interpn((Xfile, Yfile), Efile, (Xmesh, Ymesh), method='linear')
            else:
                U = interpn((Xfile, Yfile, Zfile), Ufile, (Xmesh, Ymesh, Zmesh), method='linear')
                V = interpn((Xfile, Yfile, Zfile), Vfile, (Xmesh, Ymesh, Zmesh), method='linear')
                W = interpn((Xfile, Yfile, Zfile), Wfile, (Xmesh, Ymesh, Zmesh), method='linear')
                R = interpn((Xfile, Yfile, Zfile), Rfile, (Xmesh, Ymesh, Zmesh), method='linear')
                E = interpn((Xfile, Yfile, Zfile), Efile, (Xmesh, Ymesh, Zmesh), method='linear')

            if 'changeMach' in initial_flow and initial_flow['changeMach']:
                if mesh_file['flowParameters']['Ma'] != Ma:
                    print(f"Initial flow file has a Mach number of {mesh_file['flowParameters']['Ma']} which will be changed to {Ma} to match the current simulation")
                    E *= mesh_file['flowParameters']['Ma']**2 / Ma**2

        else:
            U = flow_file['U']
            if nx != U.shape[0] or ny != U.shape[1] or (nz != U.shape[2] and U.shape[2] != 1):
                raise ValueError(f"Mesh size in initial flow file is not consistent with current mesh ({U.shape[0]}x{U.shape[1]}x{U.shape[2]} -> {nx}x{ny}x{nz})\nConsider changing the parameters file or providing a mesh file")
            V = flow_file['V']
            W = flow_file['W']
            R = flow_file['R']
            E = flow_file['E']

            if nz == 1 and U.shape[2] == 1:
                U = np.tile(U, (1, 1, nz))
                V = np.tile(V, (1, 1, nz))
                W = np.tile(W, (1, 1, nz))
                R = np.tile(R, (1, 1, nz))
                E = np.tile(E, (1, 1, nz))

    if 'addNoise' in initial_flow:
        if 'noiseVars' not in initial_flow:
            noise_vars = 'UVRE' if nz == 1 else 'UVWRE'
        else:
            noise_vars = initial_flow['noiseVars']

        noise = {}
        if 'noiseType' not in initial_flow or initial_flow['noiseType'] == 'rand':
            for var in noise_vars:
                noise[var] = initial_flow['addNoise'] * np.random.randn(nx, ny, nz)
        elif initial_flow['noiseType'] == 'uniform':
            for var in noise_vars:
                noise[var] = initial_flow['addNoise'] * np.ones((nx, ny, nz))

        if 'noiseCenter' in initial_flow:
            if 'noiseSigma' not in initial_flow:
                initial_flow['noiseSigma'] = [(X[-1] - X[0])**2 / 10, (Y[-1] - Y[0])**2 / 10, (Z[-1] - Z[0])**2 / 10]
            x0, y0, z0 = initial_flow['noiseCenter']
            sigmaX, sigmaY, sigmaZ = initial_flow['noiseSigma']
        else:
            x0, y0, z0 = 0, 0, 0
            sigmaX, sigmaY, sigmaZ = np.inf, np.inf, np.inf

        if nz == 1:
            sigmaZ = np.inf

        radius = (X[:, np.newaxis, np.newaxis] - x0)**2 / sigmaX + (Y[np.newaxis, :, np.newaxis] - y0)**2 / sigmaY
        radius += (Z[np.newaxis, np.newaxis, :] - z0)**2 / sigmaZ
        noise_gaussian = np.exp(-radius)

        for var in noise_vars:
            exec(f"{var} += noise['{var}'] * noise_gaussian")

    flow = {
        'U': U,
        'V': V,
        'W': W,
        'R': R,
        'E': E
    }

    return flow

def blasius(y):
    f0 = [0, 0, 0.33204312]
    sol = odeint(blasius_eq, f0, y)
    return sol[:, 1]

def blasius_eq(f, _):
    d = np.zeros(3)
    d[2] = -0.5 * f[0] * f[2]
    d[1] = f[2]
    d[0] = f[1]
    return d

def check_previous_run(folder):
    # Placeholder function for checking previous run
    # Implement this function based on your specific requirements
    return None



