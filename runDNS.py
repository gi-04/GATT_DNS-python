# %% Intro
# This is the main file for the DNS
# It will call all the other routines

import os
import shutil
import subprocess
import time
import numpy as np
import scipy.io as sio
from datetime import datetime

def runDNS(caseFile='parameters', extraParameters=None):
    # %% Define case file
    if caseFile is None:
        caseFile = 'parameters'

    # %% Define if simulation will actually be compiled and run or if just the preprocessing will be done
    runSimulation = True
    compileCode = True
    plotDNSDomain = False

    # %% Compiling parameters
    forceRecompileAll = False
    displayCompiling = False
    optimizeCode = True
    debugger = False
    profiler = False

    # %% Set folders for libraries
    matlabDir = ''  # Leave empty for automatic directory
    decompDir = '/usr/local/2decomp_fft'

    # %% Data logging
    logAll = False  # Save all iterations to log or just when a flow is saved

    # %% Run parameters file
    exec(open(caseFile + '.py').read())

    if extraParameters is not None:
        extraVars = unpackStruct(extraParameters)
        for var in extraVars:
            exec(f"{var} = extraParameters['{var}']")

    # %% Check if 2D or 3D
    tridimensional = mesh['z']['n'] + mesh['z']['buffer']['i']['n'] + mesh['z']['buffer']['f']['n'] > 1
    if not tridimensional:
        p_col = 1

    # %% Add source code path
    # Assuming source and source/boundaries are in the current directory
    source_path = os.path.join(os.getcwd(), 'source')
    boundaries_path = os.path.join(source_path, 'boundaries')
    os.sys.path.append(source_path)
    os.sys.path.append(boundaries_path)

    # %% Initialize log file
    if not os.path.exists(caseName):
        os.makedirs(caseName)

    print(f'Parameters file: {caseFile}\nCase name: {caseName}')

    log_file_path = os.path.join(caseName, 'log.txt')
    if not os.path.exists(log_file_path):
        with open(log_file_path, 'w') as logFile:
            logFile.write('Save number\tIteration\tSimulation time\tdt       \tCFL      \tU change\tV change\tW change\tR change\tE change')
            if 'trackedPoints' in mesh:
                for i in range(len(mesh['trackedPoints'])):
                    logFile.write(f'\tU{i}            \tV{i}            \tW{i}            \tR{i}            \tE{i}            ')
            logFile.write('\n')

    # %% Run preprocessing routine or reload previous
    bin_dir = os.path.join(caseName, 'bin')
    if not os.path.exists(bin_dir):
        os.makedirs(bin_dir)

    if forceRecompileAll:
        compileCode = True
        for file in os.listdir(bin_dir):
            if file.endswith('.mod') or file.endswith('.o'):
                os.remove(os.path.join(bin_dir, file))

    print('Running preprocessor')
    preprocessing()

    print(f'Mesh size: {mesh["nx"]} x {mesh["ny"]} x {mesh["nz"]}')

    # %% Generate initial flow if needed
    if genInitialFlow:  # If there is no previous run, generate new initial flow
        print('Generating new initial flow')
        # Compute initial flow
        flow = generateInitialFlow(mesh, flowParameters, flowType['initial'], boundary['insideWall'])

        # Save initial flow to file
        flowToSave = flow
        flowToSave['t'] = 0
        for var in 'UVWRE':
            flowToSave[var][boundary['insideWall']] = np.nan
        sio.savemat(os.path.join(caseName, 'flow_0000000000.mat'), flowToSave, do_compression=True)

        if 'meanFile' in flowType['initial']:
            flowTypeTemp = {'initial': {'type': 'file', 'flowFile': flowType['initial']['meanFile']}}
            if 'meshFile' in flowType['initial']:
                flowTypeTemp['initial']['meshFile'] = flowType['initial']['meshFile']

            meanFlow = generateInitialFlow(mesh, flowParameters, flowTypeTemp['initial'], boundary['insideWall'])

            # Save initial flow to file
            flowToSave = meanFlow
            flowToSave['t'] = 0
            for var in 'UVWRE':
                flowToSave[var][boundary['insideWall']] = np.nan
            sio.savemat(os.path.join(caseName, 'meanflowSFD.mat'), flowToSave, do_compression=True)

    else:
        print(f'Resuming from file number {time.nStep}')

    # Copy parameters file to Fortran folder and write to log2.txt
    log_file2_path = os.path.join(bin_dir, 'log2.txt')
    with open(log_file2_path, 'a') as logFile2:
        logFile2.write(f'DNS started at {datetime.now().strftime("%d-%b-%Y %H:%M:%S")}\n')
        logFile2.write(f'Parameters file: {caseFile}.py\n')
        logFile2.write(f'Starting flow file: flow_{time.nStep:010d}.mat\n\n')
        if os.path.exists(os.path.join(bin_dir, 'parameters.py')):
            parametersDiffStatus, parametersDiff = subprocess.getstatusoutput(f'diff {os.path.join(bin_dir, "parameters.py")} {caseFile}.py')
            if parametersDiffStatus:
                logFile2.write(f'Parameters file was changed:\n{parametersDiff}\n')

    shutil.copyfile(caseFile + '.py', os.path.join(bin_dir, 'parameters.py'))

    if extraParameters is not None:
        sio.savemat(os.path.join(bin_dir, 'extraParameters.mat'), extraParameters)

    # Compile code
    if compileCode:
        print('Compiling code')
        compileFortran()

    # %% Plot domain
    if plotDNSDomain:
        plotDomain()
        plt.draw()

    # %% Call Fortran code
    if runSimulation and not debugger and not profiler:
        print('Starting code')
        start_time = time.time()
        subprocess.run(f'cd {bin_dir} && mpirun -np {p_row * p_col} main {caseName}', shell=True)
        print(f"Execution time: {time.time() - start_time} seconds")
    elif runSimulation and debugger:
        print('Starting code with debugger')
        subprocess.run(f'cd {bin_dir} && mpirun -n {p_row * p_col} xterm -sl 1000000 -fg white -bg black -hold -e gdb -ex run --args ./main {caseName}', shell=True)
    elif runSimulation and profiler:
        print('Starting code with profiler')
        os.environ['GMON_OUT_PREFIX'] = 'gmon.out'
        start_time = time.time()
        subprocess.run(f'cd {bin_dir} && mpirun -np {p_row * p_col} main {caseName}', shell=True)
        print(f"Execution time: {time.time() - start_time} seconds")
        subprocess.run(f'cd {bin_dir} && gprof -l main gmon.out > profile.txt', shell=True)
        shutil.move(os.path.join(bin_dir, 'profile.txt'), '.')

    # %% Write to log2 file
    with open(log_file2_path, 'a') as logFile2:
        logFile2.write(f'DNS finished at {datetime.now().strftime("%d-%b-%Y %H:%M:%S")}\n\n')

    # %% Get outputs if needed
    flowHandles = []
    if 'flowHandles' in locals():
        allCaseFiles = os.listdir(caseName)
        for name in allCaseFiles:
            if len(name) == 19 and 'flow_' in name and name.endswith('.mat'):
                flowHandles.append(sio.loadmat(os.path.join(caseName, name)))

    info = {}
    if 'info' in locals():
        allVars = locals()
        for var in allVars:
            info[var] = allVars[var]

    return flowHandles, info

def unpackStruct(structure):
    varList = []
    for varName, value in structure.items():
        if isinstance(value, dict):
            varList2 = unpackStruct(value)
            for var in varList2:
                varList.append(f"{varName}.{var}")
        else:
            varList.append(varName)
    return varList



