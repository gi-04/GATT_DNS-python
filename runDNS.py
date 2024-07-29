# %% Intro
# This is the main file for the DNS
# It will call all the other routines

import os
import shutil
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
    # Assuming source and source/boundaries are in the Python path

    # %% Initialize log file
    if not os.path.exists(caseName):
        os.makedirs(caseName)

    print(f'Parameters file: {caseFile}\nCase name: {caseName}')

    if not os.path.exists(f'{caseName}/log.txt'):
        with open(f'{caseName}/log.txt', 'w') as logFile:
            logFile.write('Save number\tIteration\tSimulation time\tdt       \tCFL      \tU change\tV change\tW change\tR change\tE change')
            if 'trackedPoints' in mesh:
                for i in range(len(mesh['trackedPoints'])):
                    logFile.write(f'\tU{i+1}            \tV{i+1}            \tW{i+1}            \tR{i+1}            \tE{i+1}            ')
            logFile.write('\n')

    # %% Run preprocessing routine or reload previous
    if not os.path.exists(f'{caseName}/bin'):
        os.makedirs(f'{caseName}/bin')

    if forceRecompileAll:
        compileCode = True
        for file in os.listdir(f'{caseName}/bin'):
            if file.endswith('.mod') or file.endswith('.o'):
                os.remove(os.path.join(f'{caseName}/bin', file))

    print('Running preprocessor')
    preprocessing()

    print(f'Mesh size: {mesh["nx"]} x {mesh["ny"]} x {mesh["nz"]}')

    # %% Generate initial flow if needed
    if genInitialFlow:  # If there is no previous run, generate new initial flow
        print('Generating new initial flow')
        # Compute initial flow
        flow = generateInitialFlow(mesh, flowParameters, flowType['initial'], boundary['insideWall'], flowType['name'])

        # Save initial flow to file
        flowToSave = flow
        flowToSave['t'] = 0
        for var in 'UVWRE':
            flowToSave[var][boundary['insideWall']] = np.nan
        sio.savemat(f'{caseName}/flow_0000000000.mat', flowToSave, do_compression=True)

        if 'meanFile' in flowType['initial']:
            flowTypeTemp = {'initial': {'type': 'file', 'flowFile': flowType['initial']['meanFile']}}
            if 'meshFile' in flowType['initial']:
                flowTypeTemp['initial']['meshFile'] = flowType['initial']['meshFile']

            meanFlow = generateInitialFlow(mesh, flowParameters, flowTypeTemp['initial'], boundary['insideWall'], flowType['name'])

            # Save initial flow to file
            flowToSave = meanFlow
            flowToSave['t'] = 0
            for var in 'UVWRE':
                flowToSave[var][boundary['insideWall']] = np.nan
            sio.savemat(f'{caseName}/meanflowSFD.mat', flowToSave, do_compression=True)

    else:
        print(f'Resuming from file number {time["nStep"]}')

    # Copy parameters file to Fortran folder and write to log2.txt
    with open(f'{caseName}/bin/log2.txt', 'a') as logFile2:
        logFile2.write(f'DNS started at {datetime.now().strftime("%d-%b-%Y %H:%M:%S")}\n')
        logFile2.write(f'Parameters file: {caseFile}.py\n')
        logFile2.write(f'Starting flow file: flow_{time["nStep"]:010d}.mat\n\n')
        if os.path.exists(f'{caseName}/bin/parameters.py'):
            parametersDiffStatus, parametersDiff = os.system(f'diff {caseName}/bin/parameters.py {caseFile}.py')
            if parametersDiffStatus:
                logFile2.write(f'Parameters file was changed:\n{parametersDiff}\n')

    shutil.copy(f'{caseFile}.py', f'{caseName}/bin/parameters.py')

    if extraParameters is not None:
        sio.savemat(f'{caseName}/bin/extraParameters.mat', extraParameters)

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
        start_time = datetime.now()
        os.system(f'cd {caseName}/bin && mpirun -np {p_row * p_col} main {caseName}')
        end_time = datetime.now()
        print(f'Duration: {end_time - start_time}')
    elif runSimulation and debugger:
        print('Starting code with debugger')
        os.system(f'cd {caseName}/bin && mpirun -n {p_row * p_col} xterm -sl 1000000 -fg white -bg black -hold -e gdb -ex run --args ./main {caseName}')
    elif runSimulation and profiler:
        print('Starting code with profiler')
        os.system('export GMON_OUT_PREFIX="gmon.out"')
        start_time = datetime.now()
        os.system(f'cd {caseName}/bin && mpirun -np {p_row * p_col} main {caseName}')
        end_time = datetime.now()
        print(f'Duration: {end_time - start_time}')
        os.system(f'cd {caseName}/bin && gprof -l main gmon.out > profile.txt')
        shutil.move(f'{caseName}/bin/profile.txt', '.')

    # %% Write to log2 file
    with open(f'{caseName}/bin/log2.txt', 'a') as logFile2:
        logFile2.write(f'DNS finished at {datetime.now().strftime("%d-%b-%Y %H:%M:%S")}\n\n')

    # %% Get outputs if needed
    flowHandles = []
    if nargout > 0:
        allCaseFiles = os.listdir(caseName)
        for name in allCaseFiles:
            if len(name) == 19 and re.match(r'flow_\d*.mat', name):
                flowHandles.append(sio.loadmat(f'{caseName}/{name}'))

    info = {}
    if nargout == 2:
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



