import os
import re
import hashlib
import subprocess
from contextlib import contextmanager

case_name = 'ReD3000-Ddelta5-Ma03-LDinf'

files_for_mean = 90  # Number of files for each mean
spacing_for_mean = 100  # Number of saved files between each mean
n_means = 20  # Number of loops
tol = 1e-12  # Tolerance for stopping

# Load base parameters file
with open(os.path.join(case_name, 'bin', 'parameters.m'), 'r') as par_file:
    parameters = []
    for line in par_file:
        if 'time.qtimes' in line or 'time.control' in line:
            exec(line)
        parameters.append(line)

md = hashlib.md5()
md.update(case_name.encode())
hash_str = md.hexdigest()
par_file_name = f'parameters_automean_{hash_str}'

with open(f'{par_file_name}.m', 'w') as f:
    for line in parameters:
        f.write(f'{line}\n')
    f.write(f'\ncase_name = \'{case_name}\';\n')
    if time.control == 'cfl':
        f.write(f'\ntime.tmax = {t + spacing_for_mean * time.qtimes};\n')
    elif time.control == 'dt':
        f.write(f'\ntime.tmax = {lastStep + spacing_for_mean * time.qtimes};\n')

@contextmanager
def cleanup_file(file_name):
    try:
        yield
    finally:
        os.remove(f'{file_name}.m')

n_mean = 1

while True:
    changes_str = subprocess.check_output(['tail', '-1', os.path.join(case_name, 'log.txt')]).decode().strip()
    changes = [float(x) for x in changes_str.split()]
    change = max(changes[5:10])

    if n_mean > n_means or change < tol:
        break

    print(f'Computing means loop: {n_mean} Current change: {change}')
    n_mean += 1

    case_files = []
    for file_info in os.listdir(case_name):
        if len(file_info) == 19 and file_info.startswith('flow_') and file_info.endswith('.mat'):
            case_files.append(file_info)

    case_files = case_files[-files_for_mean:]

    last_step = int(case_files[-1][5:-4])

    u = v = w = r = e = 0
    for file_name in case_files:
        current = scipy.io.loadmat(os.path.join(case_name, file_name))
        u += current['U']
        v += current['V']
        w += current['W']
        r += current['R']
        e += current['E']

    u /= files_for_mean
    v /= files_for_mean
    w /= files_for_mean
    r /= files_for_mean
    e /= files_for_mean

    new_file_name = os.path.join(case_name, f'flow_{last_step + 1:010d}.mat')
    scipy.io.savemat(new_file_name, {'t': t, 'U': u, 'V': v, 'W': w, 'R': r, 'E': e})

    subprocess.run(['mv', os.path.join(case_name, 'meanflowSFD.mat'), os.path.join(case_name, f'meanflowSFD_{last_step}.mat')], stderr=subprocess.DEVNULL)

    with open(f'{par_file_name}.m', 'w') as f:
        for line in parameters:
            f.write(f'{line}\n')
        f.write(f'\ncase_name = \'{case_name}\';\n')
        if time.control == 'cfl':
            f.write(f'\ntime.tmax = {t + spacing_for_mean * time.qtimes};\n')
        elif time.control == 'dt':
            f.write(f'\ntime.tmax = {last_step + spacing_for_mean * time.qtimes};\n')

    with open(os.path.join(case_name, 'bin', 'log2.txt'), 'a') as log_file:
        log_file.write(f'Mean flow calculated from {case_files[0]} to {case_files[-1]}\n')

    subprocess.run(['rehash', 'PATH'])
    run_dns(par_file_name)

with cleanup_file(par_file_name):
    pass


