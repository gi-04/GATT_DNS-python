import os
import re
import scipy.io

def checkPreviousRun(case_name):
    """
    This function checks for previous run files.
    case_name is the folder to be checked.
    nStep is the last time step found.
    nx, ny, and nz are the size of the mesh in the saved file.
    If no files are found, empty arrays are returned.
    """
    
    allFiles = os.listdir(case_name)  # List all files

    caseFiles = []  # Flow file will be placed here

    for name in allFiles:
        if len(name) == 19 and re.search(r'flow_\d*.mat', name):  # Check the file name
            caseFiles.append(name)

    if not caseFiles:
        return [], [], [], []

    nSteps = [int(re.search(r'\d+', name).group()) for name in caseFiles]

    nStep = max(nSteps)

    if nStep:
        filePath = os.path.join(case_name, f'flow_{nStep:010d}.mat')
        fileObject = scipy.io.loadmat(filePath)
        nx, ny, nz = fileObject['U'].shape
        return nStep, nx, ny, nz

    return nStep, [], [], []
