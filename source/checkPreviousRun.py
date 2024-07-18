import os
import re
import scipy.io

def check_previous_run(case_name):
    """
    This function checks for previous run files.
    case_name is the folder to be checked.
    nStep is the last time step found.
    nx, ny, and nz are the size of the mesh in the saved file.
    If no files are found, empty arrays are returned.
    """
    
    all_files = os.listdir(case_name)  # List all files

    case_files = []  # Flow file will be placed here

    for name in all_files:
        if len(name) == 19 and re.search(r'flow_\d*.mat', name):  # Check the file name
            case_files.append(name)

    if not case_files:
        return [], [], [], []

    n_steps = [int(re.search(r'\d+', name).group()) for name in case_files]

    n_step = max(n_steps)

    if n_step:
        file_path = os.path.join(case_name, f'flow_{n_step:010d}.mat')
        file_object = scipy.io.loadmat(file_path)
        nx, ny, nz = file_object['U'].shape
        return n_step, nx, ny, nz

    return n_step, [], [], []

# Example usage:
# nStep, nx, ny, nz = check_previous_run('case_folder')



