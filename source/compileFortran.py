# This script compiles the Fortran files

import os
import subprocess

# Create extra makefile with the directories that are specific to this run
out_file_path = os.path.join(caseName, 'bin', 'makefile_extra')
with open(out_file_path, 'w') as out_file:
    if not matlabDir:
        matlabDir = os.getenv('MATLABROOT', '')

    out_file.write(f'MATROOT = {matlabDir}\n')
    out_file.write(f'DECOMPDIR = {decompDir}\n')

    if optimizeCode and not debugger:  # Optimization options
        out_file.write('ARGS += -O5 -fcheck=all -fno-finite-math-only -march=native\n')

    if debugger:  # Debugging options
        out_file.write('ARGS += -O0 -g -fbounds-check\n')
    elif profiler:  # Profiling options
        out_file.write('ARGS += -g -pg\n')

# Run make
suppress_output = '' if displayCompiling else ' >/dev/null'  # This suppresses the compiler output

main_binary_path = os.path.join(caseName, 'bin', 'main')
if os.path.exists(main_binary_path):  # Remove the main binary to force recompiling
    os.remove(main_binary_path)

make_command = f'cd {os.path.join(caseName, "bin")} && make --makefile=../../source/Fortran/makefile{suppress_output}'
status = subprocess.call(make_command, shell=True)  # Call make

if status != 0:
    raise RuntimeError('Fortran compiling has failed')



