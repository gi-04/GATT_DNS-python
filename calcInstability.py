import os
import numpy as np
import scipy.io
import scipy.sparse.linalg
import subprocess
from scipy import linalg

# Define parameters
case_folder = 'Red600-D5-Ma04-LD3-base'
comment = ''

M = 1000  # Size of Krylov space
epsilon = 1e-6  # Disturbance magnitude
dT = 5
single_beta = True
first_order = True

remove_buffer_zone = False

resume_SFD = False

# SFD parameters (optional)
# SFD.type = 0  # 0 = off, 1 = whole domain, 2 = buffer zone only;
# SFD.Delta = inf
# SFD.X = 0.005

# SFD.applyY = True

# SFD.extraRegion{1}.location = [750 0 0]
# SFD.extraRegion{1}.size = [250 5 inf]
# SFD.extraRegion{1}.X = 0.02

# Domain parameters (optional)
# nz = 16
# beta = np.pi

# Domain decomposition (optional)
# global_p_row = 4
# global_p_col = 1

# Output parameters
n_modes = 50  # Number of modes to be saved
remove_CC = True
n_krilov = [1000]  # When the results files will be generated
save_every = 100  # When intermediary files will be saved for resuming

print_eigs_every = 10  # When intermediary eigenvalues will be printed
print_eigs_n = 6  # Number of eigenvalues printed

display_all_outputs = False

# Set folders for libraries
matlab_dir = ''  # Leave empty for automatic directory
decomp_dir = '/usr/local/2decomp_fft'
use_devshm = True  # Create a symlink for the dev/shm folder for temp files

# Parameters required by the DNS
log_all = False
display_compiling = False
optimize_code = True
debugger = False
profiler = False
running_LST = True

if display_all_outputs:
    suppress_output = ''
else:
    suppress_output = ' > /dev/null 2>&1'

# Read parameters and change some values
case_file = 'parameters'
print('Reading parameters file: ')
if os.path.exists(os.path.join(case_folder, case_file + '.m')):  # First check case dir
    print(os.path.join(case_folder, case_file + '.m'))
    os.chdir(case_folder)
    exec(open(case_file).read())
    os.chdir('..')
elif os.path.exists(os.path.join(case_folder, 'bin', case_file + '.m')):  # Then check bin dir
    print(os.path.join(case_folder, 'bin', case_file + '.m'))
    os.chdir(case_folder)
    os.chdir('bin')
    exec(open(case_file).read())
    os.chdir('..')
    os.chdir('..')
else:  # Now check for any .m in case dir
    all_files = os.listdir(case_folder)
    par_files = []

    for i, name in enumerate(all_files):
        if name.endswith('.m'):
            par_files.append(name)  # #ok<AGROW>

    if len(par_files) == 1:
        print(os.path.join(case_folder, par_files[0]))
        os.chdir(case_folder)
        exec(open(par_files[0][:-2]).read())
        os.chdir('..')

    else:  # Finally, check for any .m in bin dir
        all_files = os.listdir(os.path.join(case_folder, 'bin'))
        par_files = []

        for i, name in enumerate(all_files):
            if name.endswith('.m'):
                par_files.append(name)  # #ok<AGROW>

        if len(par_files) == 1:
            print(os.path.join(case_folder, 'bin', par_files[0]))
            os.chdir(case_folder)
            os.chdir('bin')
            exec(open(par_files[0][:-2]).read())
            os.chdir('..')
            os.chdir('..')
        else:
            raise ValueError('Parameters file not found')

case_name = case_folder
log_all = False

if 'beta' in locals() and beta > 0:
    domain.zi = 0
    domain.zf = 2 * np.pi / beta
    mesh.z.fixPeriodicDomainSize = True
    if single_beta:
        mesh.z.n = 4

if 'nz' in locals():
    mesh.z.n = nz

if 'SFD' in locals():
    extra_vars = [var for var in dir() if var.startswith('SFD.')]
    for var in extra_vars:
        exec(f'numMethods.SFD.{var.split(".")[-1]} = SFD.{var};')

if 'global_p_row' in locals():
    p_row = global_p_row
if 'global_p_col' in locals():
    p_col = global_p_col

time.control = 'cfl'
time.qtimes = dT
time.tmax = dT

# Create case name
if mesh.z.n == 1:
    single_beta = False

case_name_instability = f"{case_folder}-dT{dT}-eps{epsilon:.2e}"
if not first_order:
    case_name_instability += '-SO'
if mesh.z.n > 1 and not single_beta:
    case_name_instability += '-MB'
if 'nz' in locals():
    case_name_instability += f'-nz{nz}'
if 'beta' in locals():
    case_name_instability += f'-beta{beta/np.pi:.2f}pi'
if 'comment' in locals() and comment:
    case_name_instability += f'-{comment}'
print(f'Case name: {case_name_instability}')

# Run preprocessing
print('Running preprocessor')
import source
import source.boundaries

# Set maximum number of threads for Matlab processes
os.environ['MKL_NUM_THREADS'] = str(p_row * p_col)

# Clear previous system calls cache
subprocess.run(['', ''], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# Create a copy of the original bin folder and mesh file to be restored after compling
cleanup_obj = lambda: clean_up_function(case_folder)
subprocess.run(f"cp -r {case_folder}/bin {case_folder}/bin_backup >/dev/null 2>&1", shell=True)
subprocess.run(f"cp {case_folder}/mesh.mat {case_folder}/mesh.mat.backup >/dev/null 2>&1", shell=True)

preprocessing()

# Fix nSave in parameters.F90 and set to 0
subprocess.run(f"sed -i 's/    integer :: nSave = .*/    integer :: nSave = 0/' {case_folder}/bin/parameters.F90", shell=True)

# Fix resumeSFD in parameters.F90
if resume_SFD:
    subprocess.run(f"sed -i 's/    integer :: resumeMeanFlow = .*/    integer :: resumeMeanFlow = 1/' {case_folder}/bin/parameters.F90", shell=True)

# Compile code
print('Compiling code')
subprocess.run(f"cd {case_folder}/bin && make {case_folder}", shell=True)

# Create new folders and copy files
print('Creating new folders and files')

if not os.path.exists(os.path.join(case_folder, 'Instability')):
    if use_devshm:
        inst_folder_name = f"/dev/shm/Instability{np.random.randint(1e8)}"
        os.makedirs(inst_folder_name)
        subprocess.run(f"ln -s {inst_folder_name} {case_folder}/Instability", shell=True)
    else:
        os.makedirs(os.path.join(case_folder, 'Instability'))

if not os.path.exists(os.path.join(case_folder, 'Instability', case_name_instability)):
    os.makedirs(os.path.join(case_folder, 'Instability', case_name_instability))
else:
    subprocess.run(f"rm {case_folder}/Instability/{case_name_instability}/* -r", shell=True)

subprocess.run(f"cp -r {case_folder}/bin {case_folder}/Instability/{case_name_instability}/bin", shell=True)

log_file = open(os.path.join(case_folder, 'Instability', case_name_instability, 'log.txt'), 'w')
log_file.close()

# Restore the original files
cleanup_obj = lambda: clean_up_function(case_folder, case_name_instability)

# Read the base flow
print('Reading base flow: ')

if os.path.exists(os.path.join(case_folder, 'baseflow.mat')):
    baseflow = scipy.io.loadmat(os.path.join(case_folder, 'baseflow.mat'))
    print(f"{case_folder}/baseflow.mat")
else:
    n_step = check_previous_run(case_folder)
    baseflow = scipy.io.loadmat(f"{case_folder}/flow_{n_step:010d}.mat")
    print(f"{case_folder}/flow_{n_step:010d}.mat")

# If needed, replicate flow to make it 3D
if mesh.z.n > 1 and baseflow['U'].shape[2] == 1:
    baseflow['U'] = np.repeat(baseflow['U'][:, :, np.newaxis], mesh.z.n, axis=2)
    baseflow['V'] = np.repeat(baseflow['V'][:, :, np.newaxis], mesh.z.n, axis=2)
    baseflow['W'] = np.repeat(baseflow['W'][:, :, np.newaxis], mesh.z.n, axis=2)
    baseflow['R'] = np.repeat(baseflow['R'][:, :, np.newaxis], mesh.z.n, axis=2)
    baseflow['E'] = np.repeat(baseflow['E'][:, :, np.newaxis], mesh.z.n, axis=2)

# Get the physical flow region
flow_region = ~np.isnan(baseflow['U'])  # Get region outside of walls

if remove_buffer_zone:
    flow_region[[0:mesh.x.buffer.i.n, -mesh.x.buffer.f.n:-1], :, :] = False  # Remove buffer zones
    flow_region[:, [0:mesh.y.buffer.i.n, -mesh.y.buffer.f.n:-1], :] = False
    flow_region[:, :, [0:mesh.z.buffer.i.n, -mesh.z.buffer.f.n:-1]] = False

if mesh.z.n == 1:
    N = 4 * np.sum(flow_region)
elif not single_beta:
    N = 5 * np.sum(flow_region)
else:
    N = 5 * np.sum(np.sum(flow_region[:, :, 0], axis=(0, 1)))

# Allocate variables
zeta = np.zeros((N, M))
H = np.zeros((M, M))

# Check for previous results and reload data
if os.path.exists(os.path.join(case_folder, f'results-{case_name_instability}.mat')):
    print('Previous results file found, resuming')
    previous_results = scipy.io.loadmat(os.path.join(case_folder, f'results-{case_name_instability}.mat'))

    M0 = min(previous_results['M'], M)
    zeta[:, :M0] = previous_results['zeta'][:, :M0]
    H[0:M0, 0:M0 - 1] = previous_results['H'][0:M0, 0:M0 - 1]
else:
    M0 = 1

    zeta[:, 0] = calc_initial_disturb(mesh.X, mesh.Y, mesh.Z, single_beta, flow_region)

# Set the initial SFD state
if resume_SFD:
    t = 0
    U = baseflow['U']
    V = baseflow['V']
    W = baseflow['W']
    R = baseflow['R']
    E = baseflow['E']
    scipy.io.savemat(os.path.join(case_folder, 'Instability', case_name_instability, 'meanflowSFD.mat'), {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t}, oned_as='column', do_compression=True)

# Iterate
if first_order:  # If using a first order approximation, compute the baseflow drift
    print('Iteration 0')
    t = 0

    U = baseflow['U']
    V = baseflow['V']
    W = baseflow['W']
    R = baseflow['R']
    E = baseflow['E']
    scipy.io.savemat(os.path.join(case_folder, 'Instability', case_name_instability, 'flow_0000000000.mat'), {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t}, oned_as='column', do_compression=True)

    subprocess.run(f"(cd {case_folder}/Instability/{case_name_instability}/bin && mpirun -np {p_row*p_col} main {case_name_instability}) {suppress_output}", shell=True)
    flow_minus = scipy.io.loadmat(os.path.join(case_folder, 'Instability', case_name_instability, 'flow_0000000001.mat'))
    os.remove(os.path.join(case_folder, 'Instability', case_name_instability, 'flow_0000000001.mat'))

    UM = flow2vec(flow_minus, flow_region, single_beta)

tic = time.time()
elapsed_time = 0
for k in range(M0, M + 1):
    print()
    print(f"Iteration {k} of {M}")

    # Get current disturbance
    disturbance = vec2flow(zeta[:, k - 1] * epsilon * np.sqrt(N), flow_region, single_beta)

    # Save flows, run DNS and read results
    t = 0

    U = baseflow['U'] + disturbance['U']
    V = baseflow['V'] + disturbance['V']
    W = baseflow['W'] + disturbance['W']
    R = baseflow['R'] + disturbance['R']
    E = baseflow['E'] + disturbance['E']
    scipy.io.savemat(os.path.join(case_folder, 'Instability', case_name_instability, 'flow_0000000000.mat'), {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t}, oned_as='column', do_compression=True)

    subprocess.run(f"(cd {case_folder}/Instability/{case_name_instability}/bin && mpirun -np {p_row*p_col} main {case_name_instability}) {suppress_output}", shell=True)
    flow_plus = scipy.io.loadmat(os.path.join(case_folder, 'Instability', case_name_instability, 'flow_0000000001.mat'))
    os.remove(os.path.join(case_folder, 'Instability', case_name_instability, 'flow_0000000001.mat'))

    UP = flow2vec(flow_plus, flow_region, single_beta)

    if not first_order:
        U = baseflow['U'] - disturbance['U']
        V = baseflow['V'] - disturbance['V']
        W = baseflow['W'] - disturbance['W']
        R = baseflow['R'] - disturbance['R']
        E = baseflow['E'] - disturbance['E']
        scipy.io.savemat(os.path.join(case_folder, 'Instability', case_name_instability, 'flow_0000000000.mat'), {'U': U, 'V': V, 'W': W, 'R': R, 'E': E, 't': t}, oned_as='column', do_compression=True)

        subprocess.run(f"(cd {case_folder}/Instability/{case_name_instability}/bin && mpirun -np {p_row*p_col} main {case_name_instability}) {suppress_output}", shell=True)
        flow_minus = scipy.io.loadmat(os.path.join(case_folder, 'Instability', case_name_instability, 'flow_0000000001.mat'))
        os.remove(os.path.join(case_folder, 'Instability', case_name_instability, 'flow_0000000001.mat'))

        UM = flow2vec(flow_minus, flow_region, single_beta)

        # Compute matrices
        B = (UP - UM) / (2 * epsilon * np.sqrt(N))
    else:
        B = (UP - UM) / (epsilon * np.sqrt(N))

    uPrime = B

    for j in range(1, k):
        H[j, k - 1] = np.dot(zeta[:, j - 1], B)
        uPrime -= H[j, k - 1] * zeta[:, j - 1]

    if k < M:
        H[k, k - 1] = np.linalg.norm(uPrime)
        zeta[:, k] = uPrime / H[k, k - 1]

    if np.any(np.isnan(H)) or np.any(np.isinf(H)):
        raise ValueError('DNS solution has failed')

    # Save results to file
    if k % save_every == 0 and k != M0:
        actual_M = M
        actual_H = H
        actual_zeta = zeta

        M = k

        H[M:, :] = []
        H[:, M:] = []
        zeta[:, M:] = []

        print(f"Saving results to results-{case_name_instability}.mat")
        scipy.io.savemat(os.path.join(case_folder, f'results-{case_name_instability}.mat'), {'zeta': zeta, 'H': H, 'M': M}, oned_as='column', do_compression=True)

        zeta = actual_zeta
        H = actual_H
        M = actual_M

    # Compute modes
    if k in n_krilov:
        print('Computing modes')

        [psiH, lambdaH] = scipy.linalg.eigh(H[:k, :k])
        lambdaB = np.diag(lambdaH)

        lambdaA = np.log(lambdaB) / dT
        [~, index] = np.argsort(np.real(lambdaA), axis=0)[::-1]
        lambdaA = lambdaA[index]

        psiH = psiH[:, index]

        psiA = np.dot(zeta[:, :k], psiH)

        modes = vec2modes(psiA, flow_region, single_beta, min(n_modes, k), lambdaA, remove_CC)
        modes['lambda'] = lambdaA
        modes['X'] = mesh.X
        modes['Y'] = mesh.Y
        modes['Z'] = mesh.Z

        # Save modes
        print(f"Saving modes to modes-{case_name_instability}-M{k}.mat")
        scipy.io.savemat(os.path.join(case_folder, f'modes-{case_name_instability}-M{k}.mat'), {'modes': modes}, oned_as='column', do_compression=True)

    # Print remaining time
    elapsed_time = time.time() - tic
    remaining_time = elapsed_time / k * (M - k)
    print(f"Elapsed time: {elapsed_time:.2f} s")
    print(f"Estimated time remaining: {remaining_time:.2f} s")

    # Print leading eigenvalues
    if k % print_eigs_every == 0:
        leading_eigs = scipy.sparse.linalg.eigs(H, k=min(k, print_eigs_n), which='LM', tol=1e-4, maxiter=max(300, 2 * H.shape[0]), sigma=0)
        leading_eigs = np.log(leading_eigs) / dT
        print(f"Leading eigenvalues : {leading_eigs}")
        if 'leading_eigs_prev' in locals():
            print(f"Change : {np.abs((leading_eigs[:len(leading_eigs_prev)] - leading_eigs_prev) / leading_eigs[:len(leading_eigs_prev)]).tolist()}")
        leading_eigs_prev = leading_eigs

# Extra functions

def flow2vec(flow, flow_region, single_beta):
    three_d = flow_region.shape[2] > 1
    
    if not three_d:  # 2D flow
        vec = np.concatenate([
            flow.U[flow_region],
            flow.V[flow_region],
            flow.R[flow_region],
            flow.E[flow_region]
        ])
    elif not single_beta:  # Fully 3D flow
        vec = np.concatenate([
            flow.U[flow_region],
            flow.V[flow_region],
            flow.W[flow_region],
            flow.R[flow_region],
            flow.E[flow_region]
        ])
    else:  # Single beta, 3D flow
        if flow.U.shape[2] > 1:
            nz = flow_region.shape[2]
            z_temp = np.linspace(0, 2*np.pi, nz+1)[:-1]
            z_temp = z_temp.reshape(1, 1, -1)

            cos_z = np.cos(z_temp) * 2 / nz
            sin_z = np.sin(z_temp) * 2 / nz
            
            flow.U = np.sum(cos_z * flow.U, axis=2)
            flow.V = np.sum(cos_z * flow.V, axis=2)
            flow.W = np.sum(sin_z * flow.W, axis=2)
            flow.R = np.sum(cos_z * flow.R, axis=2)
            flow.E = np.sum(cos_z * flow.E, axis=2)
        
        flow_region = flow_region[:,:,0]
        
        vec = np.concatenate([
            flow.U[flow_region],
            flow.V[flow_region],
            flow.W[flow_region],
            flow.R[flow_region],
            flow.E[flow_region]
        ])
    
    return vec

def vec2flow(vec, flow_region, single_beta):
    nx, ny, nz = flow_region.shape
    N = np.sum(flow_region)
    
    three_d = nz > 1
    
    flow = type('Flow', (), {})()
    
    if not three_d:  # 2D flow
        flow.U = np.zeros((nx, ny))
        flow.V = np.zeros((nx, ny))
        flow.W = np.zeros((nx, ny))
        flow.R = np.zeros((nx, ny))
        flow.E = np.zeros((nx, ny))
        
        flow.U[flow_region] = vec[:N]
        flow.V[flow_region] = vec[N:2*N]
        flow.R[flow_region] = vec[2*N:3*N]
        flow.E[flow_region] = vec[3*N:4*N]
        
    elif not single_beta:  # Fully 3D flow
        flow.U = np.zeros((nx, ny, nz))
        flow.V = np.zeros((nx, ny, nz))
        flow.W = np.zeros((nx, ny, nz))
        flow.R = np.zeros((nx, ny, nz))
        flow.E = np.zeros((nx, ny, nz))
        
        flow.U[flow_region] = vec[:N]
        flow.V[flow_region] = vec[N:2*N]
        flow.W[flow_region] = vec[2*N:3*N]
        flow.R[flow_region] = vec[3*N:4*N]
        flow.E[flow_region] = vec[4*N:5*N]
    
    else:  # Single beta, 3D flow
        flow.U = np.zeros((nx, ny))
        flow.V = np.zeros((nx, ny))
        flow.W = np.zeros((nx, ny))
        flow.R = np.zeros((nx, ny))
        flow.E = np.zeros((nx, ny))
        
        N = round(N / nz)
        
        flow.U[flow_region[:,:,0]] = vec[:N]
        flow.V[flow_region[:,:,0]] = vec[N:2*N]
        flow.W[flow_region[:,:,0]] = vec[2*N:3*N]
        flow.R[flow_region[:,:,0]] = vec[3*N:4*N]
        flow.E[flow_region[:,:,0]] = vec[4*N:5*N]
    
        z_temp = np.linspace(0, 2*np.pi, nz+1)[:-1]
        z_temp = z_temp.reshape(1, 1, -1)
        
        cos_z = np.cos(z_temp)
        sin_z = np.sin(z_temp)
        
        flow.U = flow.U[:,:,np.newaxis] * cos_z
        flow.V = flow.V[:,:,np.newaxis] * cos_z
        flow.W = flow.W[:,:,np.newaxis] * sin_z
        flow.R = flow.R[:,:,np.newaxis] * cos_z
        flow.E = flow.E[:,:,np.newaxis] * cos_z
    
    return flow

def vec2modes(vec, flow_region, single_beta, n_modes, lambda_vals, remove_cc):
    nx, ny, nz = flow_region.shape
    N = np.sum(flow_region)
    
    three_d = nz > 1
    
    modes = type('Modes', (), {})()
    
    if not three_d:  # 2D flow
        modes.U = np.full((nx, ny, 1, n_modes), np.nan)
        modes.V = np.full((nx, ny, 1, n_modes), np.nan)
        modes.W = np.full((nx, ny, 1, n_modes), np.nan)
        modes.R = np.full((nx, ny, 1, n_modes), np.nan)
        modes.E = np.full((nx, ny, 1, n_modes), np.nan)
        
        U = np.full((nx, ny), np.nan)
        V = np.full((nx, ny), np.nan)
        R = np.full((nx, ny), np.nan)
        E = np.full((nx, ny), np.nan)
        
        for i in range(n_modes):
            if np.imag(lambda_vals[i]) >= 0:
                vec[:, i] = vec[:, i] / np.max(np.abs(vec[:, i]))
                
                U[flow_region] = vec[:N, i]
                V[flow_region] = vec[N:2*N, i]
                R[flow_region] = vec[2*N:3*N, i]
                E[flow_region] = vec[3*N:4*N, i]
                
                modes.U[:,:,0,i] = U
                modes.V[:,:,0,i] = V
                modes.R[:,:,0,i] = R
                modes.E[:,:,0,i] = E
        
    elif not single_beta:  # Fully 3D flow
        modes.U = np.full((nx, ny, nz, n_modes), np.nan)
        modes.V = np.full((nx, ny, nz, n_modes), np.nan)
        modes.W = np.full((nx, ny, nz, n_modes), np.nan)
        modes.R = np.full((nx, ny, nz, n_modes), np.nan)
        modes.E = np.full((nx, ny, nz, n_modes), np.nan)
        
        U = np.full((nx, ny, nz), np.nan)
        V = np.full((nx, ny, nz), np.nan)
        W = np.full((nx, ny, nz), np.nan)
        R = np.full((nx, ny, nz), np.nan)
        E = np.full((nx, ny, nz), np.nan)
        
        for i in range(n_modes):
            if np.imag(lambda_vals[i]) >= 0:
                vec[:, i] = vec[:, i] / np.max(np.abs(vec[:, i]))
                
                U[flow_region] = vec[:N, i]
                V[flow_region] = vec[N:2*N, i]
                W[flow_region] = vec[2*N:3*N, i]
                R[flow_region] = vec[3*N:4*N, i]
                E[flow_region] = vec[4*N:5*N, i]
                
                modes.U[:,:,:,i] = U
                modes.V[:,:,:,i] = V
                modes.W[:,:,:,i] = W
                modes.R[:,:,:,i] = R
                modes.E[:,:,:,i] = E
    
    else:  # Single beta, 3D flow
        N = round(N / nz)
        flow_region = flow_region[:,:,0]
        
        modes.U = np.full((nx, ny, 1, n_modes), np.nan)
        modes.V = np.full((nx, ny, 1, n_modes), np.nan)
        modes.W = np.full((nx, ny, 1, n_modes), np.nan)
        modes.R = np.full((nx, ny, 1, n_modes), np.nan)
        modes.E = np.full((nx, ny, 1, n_modes), np.nan)
        
        U = np.full((nx, ny), np.nan)
        V = np.full((nx, ny), np.nan)
        W = np.full((nx, ny), np.nan)
        R = np.full((nx, ny), np.nan)
        E = np.full((nx, ny), np.nan)
        
        for i in range(n_modes):
            if np.imag(lambda_vals[i]) >= 0:
                vec[:, i] = vec[:, i] / np.max(np.abs(vec[:, i]))
                
                U[flow_region] = vec[:N, i]
                V[flow_region] = vec[N:2*N, i]
                W[flow_region] = vec[2*N:3*N, i]
                R[flow_region] = vec[3*N:4*N, i]
                E[flow_region] = vec[4*N:5*N, i]
                
                modes.U[:,:,0,i] = U
                modes.V[:,:,0,i] = V
                modes.W[:,:,0,i] = W
                modes.R[:,:,0,i] = R
                modes.E[:,:,0,i] = E
    
    return modes

def calc_initial_disturb(X, Y, Z, single_beta, flow_region):
    x0i = np.argmin(np.gradient(X))
    y0i = np.argmin(np.gradient(Y))
    z0i = np.argmin(np.gradient(Z))

    x0 = X[x0i]
    y0 = Y[y0i]
    z0 = Z[z0i]
    
    alpha_x = 10 / (X[-1] - X[0])**2
    alpha_y = 10 / (Y[-1] - Y[0])**2
    alpha_z = 10 / (Z[-1] - Z[0])**2
    
    flow = type('Flow', (), {})()
    
    if len(Z) == 1:
        eta = np.exp(-(alpha_x*(X[:,np.newaxis]-x0)**2 + alpha_y*(Y[np.newaxis,:]-y0)**2))
        
        flow.U = eta
        flow.V = eta
        flow.W = eta
        flow.R = eta
        flow.E = eta
        
        zeta = flow2vec(flow, flow_region, single_beta)
        
    elif not single_beta:
        eta = np.exp(-(alpha_x*(X[:,np.newaxis,np.newaxis]-x0)**2 + 
                       alpha_y*(Y[np.newaxis,:,np.newaxis]-y0)**2 + 
                       alpha_z*(Z[np.newaxis,np.newaxis,:]-z0)**2))
        
        flow.U = eta
        flow.V = eta
        flow.W = eta
        flow.R = eta
        flow.E = eta
        
        zeta = flow2vec(flow, flow_region, single_beta)
        
    else:
        eta = np.exp(-(alpha_x*(X[:,np.newaxis]-x0)**2 + alpha_y*(Y[np.newaxis,:]-y0)**2))
        
        flow.U = eta
        flow.V = eta
        flow.W = eta
        flow.R = eta
        flow.E = eta
        
        zeta = flow2vec(flow, flow_region, single_beta)
    
    zeta = zeta / linalg.norm(zeta)
    
    return zeta

def sec2str(sec):
    sec = round(sec)
    if sec < 60:
        return f"{sec}s"
    else:
        min = sec // 60
        sec = sec % 60
    if min < 60:
        return f"{min}min {sec:02d}s"
    else:
        hour = min // 60
        min = min % 60
    if hour < 24:
        return f"{hour}h {min:02d}min {sec:02d}s"
    else:
        day = hour // 24
        hour = hour % 24
    return f"{day}d {hour:02d}h {min:02d}min {sec:02d}s"

def unpack_struct(structure):
    var_list = []
    for var_name, value in structure.__dict__.items():
        if isinstance(value, type):
            var_list2 = unpack_struct(value)
            var_list.extend([f"{var_name}.{item}" for item in var_list2])
        else:
            var_list.append(var_name)
    return var_list

import os
import shutil

def clean_up_function(case_folder, case_name_instability=None):
    if case_name_instability is None:  # Restore the original files
        shutil.move(os.path.join(case_folder, 'mesh.mat.backup'), 
                    os.path.join(case_folder, 'mesh.mat'))
        shutil.rmtree(os.path.join(case_folder, 'bin'), ignore_errors=True)
        shutil.move(os.path.join(case_folder, 'bin_backup'), 
                    os.path.join(case_folder, 'bin'))
    else:  # Clean up instability folder
        shutil.rmtree(os.path.join(case_folder, 'Instability', case_name_instability))
        instability_dir = os.path.join(case_folder, 'Instability')
        if len(os.listdir(instability_dir)) == 0:  # Check if Instability dir is empty and erase it
            os.remove(os.readlink(instability_dir))
            os.rmdir(instability_dir)


