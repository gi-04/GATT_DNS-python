import numpy as np

def write_fortran_matrices(case_name, matrices, num_methods, mesh):
    """
    This function writes matrices.F90, which contains all the matrices for derivatives and filters,
    as well as the regions to which they will be applied. The coefficients for Neumann boundary
    conditions are also defined here.
    """
    out_file = open(f"{case_name}/bin/matrices.F90", 'w')

    I = matrices['y']['types'].shape[0]
    J = matrices['x']['types'].shape[0]
    K = matrices['x']['types'].shape[1]
    n_procs = len(matrices['x']['blocks'])

    # Write blocks
    out_file.write('    select case(nrank)\n')

    for i in range(n_procs):
        out_file.write(f'        case({i})\n')

        # For X
        blocks = matrices['x']['blocks'][i]
        n_blocks = blocks.shape[0]
        out_file.write(f'            nDerivBlocksX = {n_blocks}\n')
        out_file.write(f'            allocate(derivBlocksX({n_blocks},5))\n')

        out_file.write('            derivBlocksX = reshape((/')
        out_file.write(','.join(map(str, blocks.flatten())))
        out_file.write('/),shape(derivBlocksX))\n\n')

        # For Y
        blocks = matrices['y']['blocks'][i]
        n_blocks = blocks.shape[0]
        out_file.write(f'            nDerivBlocksY = {n_blocks}\n')
        out_file.write(f'            allocate(derivBlocksY({n_blocks},5))\n')

        out_file.write('            derivBlocksY = reshape((/')
        out_file.write(','.join(map(str, blocks.flatten())))
        out_file.write('/),shape(derivBlocksY))\n\n')

        # For Z
        if K > 1:
            blocks = matrices['z']['blocks'][i]
            n_blocks = blocks.shape[0]
            out_file.write(f'            nDerivBlocksZ = {n_blocks}\n')
            out_file.write(f'            allocate(derivBlocksZ({n_blocks},5))\n')

            out_file.write('            derivBlocksZ = reshape((/')
            out_file.write(','.join(map(str, blocks.flatten())))
            out_file.write('/),shape(derivBlocksZ))\n\n')

    out_file.write('    end select\n\n')

    # Write filter info
    out_file.write(f'    filterX = {num_methods["filterDirections"][0]}\n')
    out_file.write(f'    filterY = {num_methods["filterDirections"][1]}\n')
    out_file.write(f'    filterZ = {num_methods["filterDirections"][2]}\n')

    # Write matrices
    # For X derivatives
    out_file.write(f'    derivnRHSx = {matrices["x"]["nRHS"]}\n')
    out_file.write(f'    filternRHSx = {matrices["x"]["nRHSf"]}\n')

    out_file.write(f'    allocate(derivsAX({I-1},{matrices["x"]["nTypes"]}))\n')
    out_file.write(f'    allocate(derivsBX({I-matrices["x"]["periodic"]},{matrices["x"]["nTypes"]}))\n')
    out_file.write(f'    allocate(derivsCX({I-1},{matrices["x"]["nTypes"]}))\n')
    out_file.write(f'    allocate(derivsRX({I},{2*matrices["x"]["nRHS"]-1},{matrices["x"]["nTypes"]}))\n\n')

    out_file.write(f'    allocate(filterAX({I-1},{matrices["x"]["nTypes"]}))\n')
    out_file.write(f'    allocate(filterBX({I-matrices["x"]["periodic"]},{matrices["x"]["nTypes"]}))\n')
    out_file.write(f'    allocate(filterCX({I-1},{matrices["x"]["nTypes"]}))\n')
    out_file.write(f'    allocate(filterRX({I},{2*matrices["x"]["nRHSf"]-1},{matrices["x"]["nTypes"]}))\n\n')

    write_matrix(out_file, 'derivsAX', matrices['x']['A'])
    write_matrix(out_file, 'derivsBX', matrices['x']['B'])
    write_matrix(out_file, 'derivsCX', matrices['x']['C'])
    write_matrix(out_file, 'derivsRX', matrices['x']['R'])

    write_matrix(out_file, 'filterAX', matrices['x']['Af'])
    write_matrix(out_file, 'filterBX', matrices['x']['Bf'])
    write_matrix(out_file, 'filterCX', matrices['x']['Cf'])
    write_matrix(out_file, 'filterRX', matrices['x']['Rf'])

    out_file.write(f'    periodicX = {matrices["x"]["periodic"]}\n\n')

    if matrices['x']['periodic']:
        out_file.write(f'    allocate(derivsDX({I},{matrices["x"]["nTypes"]}))\n')
        out_file.write(f'    allocate(filterDX({I},{matrices["x"]["nTypes"]}))\n')

        write_matrix(out_file, 'derivsDX', matrices['x']['D'])
        write_matrix(out_file, 'filterDX', matrices['x']['Df'])

    # For Y derivatives
    out_file.write(f'    derivnRHSy = {matrices["y"]["nRHS"]}\n')
    out_file.write(f'    filternRHSy = {matrices["y"]["nRHSf"]}\n')

    out_file.write(f'    allocate(derivsAY({J-1},{matrices["y"]["nTypes"]}))\n')
    out_file.write(f'    allocate(derivsBY({J-matrices["y"]["periodic"]},{matrices["y"]["nTypes"]}))\n')
    out_file.write(f'    allocate(derivsCY({J-1},{matrices["y"]["nTypes"]}))\n')
    out_file.write(f'    allocate(derivsRY({J},{2*matrices["y"]["nRHS"]-1},{matrices["y"]["nTypes"]}))\n\n')

    out_file.write(f'    allocate(filterAY({J-1},{matrices["y"]["nTypes"]}))\n')
    out_file.write(f'    allocate(filterBY({J-matrices["y"]["periodic"]},{matrices["y"]["nTypes"]}))\n')
    out_file.write(f'    allocate(filterCY({J-1},{matrices["y"]["nTypes"]}))\n')
    out_file.write(f'    allocate(filterRY({J},{2*matrices["y"]["nRHSf"]-1},{matrices["y"]["nTypes"]}))\n\n')

    write_matrix(out_file, 'derivsAY', matrices['y']['A'])
    write_matrix(out_file, 'derivsBY', matrices['y']['B'])
    write_matrix(out_file, 'derivsCY', matrices['y']['C'])
    write_matrix(out_file, 'derivsRY', matrices['y']['R'])

    write_matrix(out_file, 'filterAY', matrices['y']['Af'])
    write_matrix(out_file, 'filterBY', matrices['y']['Bf'])
    write_matrix(out_file, 'filterCY', matrices['y']['Cf'])
    write_matrix(out_file, 'filterRY', matrices['y']['Rf'])

    out_file.write(f'    periodicY = {matrices["y"]["periodic"]}\n\n')

    if matrices['y']['periodic']:
        out_file.write(f'    allocate(derivsDY({J},{matrices["y"]["nTypes"]}))\n')
        out_file.write(f'    allocate(filterDY({J},{matrices["y"]["nTypes"]}))\n')

        write_matrix(out_file, 'derivsDY', matrices['y']['D'])
        write_matrix(out_file, 'filterDY', matrices['y']['Df'])

    # For Z derivatives
    if K > 1:
        out_file.write(f'    derivnRHSz = {matrices["z"]["nRHS"]}\n')
        out_file.write(f'    filternRHSz = {matrices["z"]["nRHSf"]}\n')

        out_file.write(f'    allocate(derivsAZ({K-1},{matrices["z"]["nTypes"]}))\n')
        out_file.write(f'    allocate(derivsBZ({K-matrices["z"]["periodic"]},{matrices["z"]["nTypes"]}))\n')
        out_file.write(f'    allocate(derivsCZ({K-1},{matrices["z"]["nTypes"]}))\n')
        out_file.write(f'    allocate(derivsRZ({K},{2*matrices["z"]["nRHS"]-1},{matrices["z"]["nTypes"]}))\n\n')

        out_file.write(f'    allocate(filterAZ({K-1},{matrices["z"]["nTypes"]}))\n')
        out_file.write(f'    allocate(filterBZ({K-matrices["z"]["periodic"]},{matrices["z"]["nTypes"]}))\n')
        out_file.write(f'    allocate(filterCZ({K-1},{matrices["z"]["nTypes"]}))\n')
        out_file.write(f'    allocate(filterRZ({K},{2*matrices["z"]["nRHSf"]-1},{matrices["z"]["nTypes"]}))\n\n')

        write_matrix(out_file, 'derivsAZ', matrices['z']['A'])
        write_matrix(out_file, 'derivsBZ', matrices['z']['B'])
        write_matrix(out_file, 'derivsCZ', matrices['z']['C'])
        write_matrix(out_file, 'derivsRZ', matrices['z']['R'])

        write_matrix(out_file, 'filterAZ', matrices['z']['Af'])
        write_matrix(out_file, 'filterBZ', matrices['z']['Bf'])
        write_matrix(out_file, 'filterCZ', matrices['z']['Cf'])
        write_matrix(out_file, 'filterRZ', matrices['z']['Rf'])

        out_file.write(f'    periodicZ = {matrices["z"]["periodic"]}\n\n')

        if matrices['z']['periodic']:
            out_file.write(f'    allocate(derivsDZ({K},{matrices["z"]["nTypes"]}))\n')
            out_file.write(f'    allocate(filterDZ({K},{matrices["z"]["nTypes"]}))\n')

            write_matrix(out_file, 'derivsDZ', matrices['z']['D'])
            write_matrix(out_file, 'filterDZ', matrices['z']['Df'])

    # Neumann Coefficients
    out_file.write(f'\n    neumannLength = {len(matrices["neumannCoeffs"])}\n')
    out_file.write(f'    allocate(neumannCoeffs({len(matrices["neumannCoeffs"])}))\n')
    out_file.write('    neumannCoeffs = (/')
    out_file.write(','.join(f'{coeff:.20f}d0' for coeff in matrices['neumannCoeffs']))
    out_file.write('/)\n\n')

    out_file.write(f'\n    neumann2Length = {len(matrices["neumann2Coeffs"])}\n')
    out_file.write(f'    allocate(neumann2Coeffs({len(matrices["neumann2Coeffs"])}))\n')
    out_file.write('    neumann2Coeffs = (/')
    out_file.write(','.join(f'{coeff:.20f}d0' for coeff in matrices['neumann2Coeffs']))
    out_file.write('/)\n')

    # Tracked points
    if 'trackedPoints' not in mesh or not mesh['trackedPoints']:
        out_file.write('    nTracked = 0\n')
    else:
        out_file.write(f'    nTracked = {mesh["trackedPoints"].shape[0]}\n')
        out_file.write(f'    allocate(indTracked({mesh["trackedPoints"].shape[0]},3))\n')

        ind_tracked = mesh['trackedPoints']

        # xTemp is a vector that only contains the nodes that were originally in the mesh in case extraRefinement is used
        # This is to make sure that the tracked point remain the same when using extraRefinement
        x_temp = np.full_like(mesh['X'], np.nan)
        x_temp[::mesh['x']['extraRefinement']+1] = mesh['X'][::mesh['x']['extraRefinement']+1]
        y_temp = np.full_like(mesh['Y'], np.nan)
        y_temp[::mesh['y']['extraRefinement']+1] = mesh['Y'][::mesh['y']['extraRefinement']+1]
        z_temp = np.full_like(mesh['Z'], np.nan)
        z_temp[::mesh['z']['extraRefinement']+1] = mesh['Z'][::mesh['z']['extraRefinement']+1]

        for i in range(mesh['trackedPoints'].shape[0]):
            ind_tracked[i, 0] = np.argmin(np.abs(ind_tracked[i, 0] - x_temp))
            ind_tracked[i, 1] = np.argmin(np.abs(ind_tracked[i, 1] - y_temp))
            ind_tracked[i, 2] = np.argmin(np.abs(ind_tracked[i, 2] - z_temp))

        out_file.write('    indTracked = reshape((/')
        out_file.write(','.join(map(str, ind_tracked.flatten())))
        out_file.write('/),shape(indTracked))\n')

    out_file.close()

def write_matrix(out_file, name, var):
    out_file.write(f'    {name} = reshape((/')
    out_file.write(','.join(f'{v:.20f}d0' for v in var.flatten()))
    out_file.write(f'/),shape({name}))\n')



