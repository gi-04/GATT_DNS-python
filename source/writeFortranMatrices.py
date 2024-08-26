import numpy as np

def writeFortranMatrices(case_name, matrices, num_methods, mesh):
    """
    This function writes matrices.F90, which contains all the matrices for derivatives and filters,
    as well as the regions to which they will be applied. The coefficients for Neumann boundary
    conditions are also defined here.
    """
    outFile = open(f"{case_name}/bin/matrices.F90", 'w')

    I = matrices['y']['types'].shape[0]
    J = matrices['x']['types'].shape[0]
    K = matrices['x']['types'].shape[1]
    n_procs = len(matrices['x']['blocks'])

    # Write blocks
    outFile.write('    select case(nrank)\n')

    for i in range(n_procs):
        outFile.write(f'        case({i})\n')

        # For X
        blocks = matrices['x']['blocks'][i]
        nBlocks = blocks.shape[0]
        outFile.write(f'            nDerivBlocksX = {nBlocks}\n')
        outFile.write(f'            allocate(derivBlocksX({nBlocks},5))\n')

        outFile.write('            derivBlocksX = reshape((/')
        outFile.write(','.join(map(str, blocks.flatten())))
        outFile.write('/),shape(derivBlocksX))\n\n')

        # For Y
        blocks = matrices['y']['blocks'][i]
        nBlocks = blocks.shape[0]
        outFile.write(f'            nDerivBlocksY = {nBlocks}\n')
        outFile.write(f'            allocate(derivBlocksY({nBlocks},5))\n')

        outFile.write('            derivBlocksY = reshape((/')
        outFile.write(','.join(map(str, blocks.flatten())))
        outFile.write('/),shape(derivBlocksY))\n\n')

        # For Z
        if K > 1:
            blocks = matrices['z']['blocks'][i]
            nBlocks = blocks.shape[0]
            outFile.write(f'            nDerivBlocksZ = {nBlocks}\n')
            outFile.write(f'            allocate(derivBlocksZ({nBlocks},5))\n')

            outFile.write('            derivBlocksZ = reshape((/')
            outFile.write(','.join(map(str, blocks.flatten())))
            outFile.write('/),shape(derivBlocksZ))\n\n')

    outFile.write('    end select\n\n')

    # Write filter info
    outFile.write(f'    filterX = {num_methods["filterDirections"][0]}\n')
    outFile.write(f'    filterY = {num_methods["filterDirections"][1]}\n')
    outFile.write(f'    filterZ = {num_methods["filterDirections"][2]}\n')

    # Write matrices
    # For X derivatives
    outFile.write(f'    derivnRHSx = {matrices["x"]["nRHS"]}\n')
    outFile.write(f'    filternRHSx = {matrices["x"]["nRHSf"]}\n')

    outFile.write(f'    allocate(derivsAX({I-1},{matrices["x"]["nTypes"]}))\n')
    outFile.write(f'    allocate(derivsBX({I-matrices["x"]["periodic"]},{matrices["x"]["nTypes"]}))\n')
    outFile.write(f'    allocate(derivsCX({I-1},{matrices["x"]["nTypes"]}))\n')
    outFile.write(f'    allocate(derivsRX({I},{2*matrices["x"]["nRHS"]-1},{matrices["x"]["nTypes"]}))\n\n')

    outFile.write(f'    allocate(filterAX({I-1},{matrices["x"]["nTypes"]}))\n')
    outFile.write(f'    allocate(filterBX({I-matrices["x"]["periodic"]},{matrices["x"]["nTypes"]}))\n')
    outFile.write(f'    allocate(filterCX({I-1},{matrices["x"]["nTypes"]}))\n')
    outFile.write(f'    allocate(filterRX({I},{2*matrices["x"]["nRHSf"]-1},{matrices["x"]["nTypes"]}))\n\n')

    writeMatrix(outFile, 'derivsAX', matrices['x']['A'])
    writeMatrix(outFile, 'derivsBX', matrices['x']['B'])
    writeMatrix(outFile, 'derivsCX', matrices['x']['C'])
    writeMatrix(outFile, 'derivsRX', matrices['x']['R'])

    writeMatrix(outFile, 'filterAX', matrices['x']['Af'])
    writeMatrix(outFile, 'filterBX', matrices['x']['Bf'])
    writeMatrix(outFile, 'filterCX', matrices['x']['Cf'])
    writeMatrix(outFile, 'filterRX', matrices['x']['Rf'])

    outFile.write(f'    periodicX = {matrices["x"]["periodic"]}\n\n')

    if matrices['x']['periodic']:
        outFile.write(f'    allocate(derivsDX({I},{matrices["x"]["nTypes"]}))\n')
        outFile.write(f'    allocate(filterDX({I},{matrices["x"]["nTypes"]}))\n')

        writeMatrix(outFile, 'derivsDX', matrices['x']['D'])
        writeMatrix(outFile, 'filterDX', matrices['x']['Df'])

    # For Y derivatives
    outFile.write(f'    derivnRHSy = {matrices["y"]["nRHS"]}\n')
    outFile.write(f'    filternRHSy = {matrices["y"]["nRHSf"]}\n')

    outFile.write(f'    allocate(derivsAY({J-1},{matrices["y"]["nTypes"]}))\n')
    outFile.write(f'    allocate(derivsBY({J-matrices["y"]["periodic"]},{matrices["y"]["nTypes"]}))\n')
    outFile.write(f'    allocate(derivsCY({J-1},{matrices["y"]["nTypes"]}))\n')
    outFile.write(f'    allocate(derivsRY({J},{2*matrices["y"]["nRHS"]-1},{matrices["y"]["nTypes"]}))\n\n')

    outFile.write(f'    allocate(filterAY({J-1},{matrices["y"]["nTypes"]}))\n')
    outFile.write(f'    allocate(filterBY({J-matrices["y"]["periodic"]},{matrices["y"]["nTypes"]}))\n')
    outFile.write(f'    allocate(filterCY({J-1},{matrices["y"]["nTypes"]}))\n')
    outFile.write(f'    allocate(filterRY({J},{2*matrices["y"]["nRHSf"]-1},{matrices["y"]["nTypes"]}))\n\n')

    writeMatrix(outFile, 'derivsAY', matrices['y']['A'])
    writeMatrix(outFile, 'derivsBY', matrices['y']['B'])
    writeMatrix(outFile, 'derivsCY', matrices['y']['C'])
    writeMatrix(outFile, 'derivsRY', matrices['y']['R'])

    writeMatrix(outFile, 'filterAY', matrices['y']['Af'])
    writeMatrix(outFile, 'filterBY', matrices['y']['Bf'])
    writeMatrix(outFile, 'filterCY', matrices['y']['Cf'])
    writeMatrix(outFile, 'filterRY', matrices['y']['Rf'])

    outFile.write(f'    periodicY = {matrices["y"]["periodic"]}\n\n')

    if matrices['y']['periodic']:
        outFile.write(f'    allocate(derivsDY({J},{matrices["y"]["nTypes"]}))\n')
        outFile.write(f'    allocate(filterDY({J},{matrices["y"]["nTypes"]}))\n')

        writeMatrix(outFile, 'derivsDY', matrices['y']['D'])
        writeMatrix(outFile, 'filterDY', matrices['y']['Df'])

    # For Z derivatives
    if K > 1:
        outFile.write(f'    derivnRHSz = {matrices["z"]["nRHS"]}\n')
        outFile.write(f'    filternRHSz = {matrices["z"]["nRHSf"]}\n')

        outFile.write(f'    allocate(derivsAZ({K-1},{matrices["z"]["nTypes"]}))\n')
        outFile.write(f'    allocate(derivsBZ({K-matrices["z"]["periodic"]},{matrices["z"]["nTypes"]}))\n')
        outFile.write(f'    allocate(derivsCZ({K-1},{matrices["z"]["nTypes"]}))\n')
        outFile.write(f'    allocate(derivsRZ({K},{2*matrices["z"]["nRHS"]-1},{matrices["z"]["nTypes"]}))\n\n')

        outFile.write(f'    allocate(filterAZ({K-1},{matrices["z"]["nTypes"]}))\n')
        outFile.write(f'    allocate(filterBZ({K-matrices["z"]["periodic"]},{matrices["z"]["nTypes"]}))\n')
        outFile.write(f'    allocate(filterCZ({K-1},{matrices["z"]["nTypes"]}))\n')
        outFile.write(f'    allocate(filterRZ({K},{2*matrices["z"]["nRHSf"]-1},{matrices["z"]["nTypes"]}))\n\n')

        writeMatrix(outFile, 'derivsAZ', matrices['z']['A'])
        writeMatrix(outFile, 'derivsBZ', matrices['z']['B'])
        writeMatrix(outFile, 'derivsCZ', matrices['z']['C'])
        writeMatrix(outFile, 'derivsRZ', matrices['z']['R'])

        writeMatrix(outFile, 'filterAZ', matrices['z']['Af'])
        writeMatrix(outFile, 'filterBZ', matrices['z']['Bf'])
        writeMatrix(outFile, 'filterCZ', matrices['z']['Cf'])
        writeMatrix(outFile, 'filterRZ', matrices['z']['Rf'])

        outFile.write(f'    periodicZ = {matrices["z"]["periodic"]}\n\n')

        if matrices['z']['periodic']:
            outFile.write(f'    allocate(derivsDZ({K},{matrices["z"]["nTypes"]}))\n')
            outFile.write(f'    allocate(filterDZ({K},{matrices["z"]["nTypes"]}))\n')

            writeMatrix(outFile, 'derivsDZ', matrices['z']['D'])
            writeMatrix(outFile, 'filterDZ', matrices['z']['Df'])

    # Neumann Coefficients
    outFile.write(f'\n    neumannLength = {len(matrices["neumannCoeffs"])}\n')
    outFile.write(f'    allocate(neumannCoeffs({len(matrices["neumannCoeffs"])}))\n')
    outFile.write('    neumannCoeffs = (/')
    outFile.write(','.join(f'{coeff:.20f}d0' for coeff in matrices['neumannCoeffs']))
    outFile.write('/)\n\n')

    outFile.write(f'\n    neumann2Length = {len(matrices["neumann2Coeffs"])}\n')
    outFile.write(f'    allocate(neumann2Coeffs({len(matrices["neumann2Coeffs"])}))\n')
    outFile.write('    neumann2Coeffs = (/')
    outFile.write(','.join(f'{coeff:.20f}d0' for coeff in matrices['neumann2Coeffs']))
    outFile.write('/)\n')

    # Tracked points
    if 'trackedPoints' not in mesh or not mesh['trackedPoints']:
        outFile.write('    nTracked = 0\n')
    else:
        outFile.write(f'    nTracked = {mesh["trackedPoints"].shape[0]}\n')
        outFile.write(f'    allocate(indTracked({mesh["trackedPoints"].shape[0]},3))\n')

        ind_tracked = mesh['trackedPoints']

        # xTemp is a vector that only contains the nodes that were originally in the mesh in case extraRefinement is used
        # This is to make sure that the tracked point remain the same when using extraRefinement
        xTemp = np.full_like(mesh['X'], np.nan)
        xTemp[::mesh['x']['extraRefinement']+1] = mesh['X'][::mesh['x']['extraRefinement']+1]
        YTemp = np.full_like(mesh['Y'], np.nan)
        YTemp[::mesh['y']['extraRefinement']+1] = mesh['Y'][::mesh['y']['extraRefinement']+1]
        ZTemp = np.full_like(mesh['Z'], np.nan)
        ZTemp[::mesh['z']['extraRefinement']+1] = mesh['Z'][::mesh['z']['extraRefinement']+1]

        for i in range(mesh['trackedPoints'].shape[0]):
            ind_tracked[i, 0] = np.argmin(np.abs(ind_tracked[i, 0] - xTemp))
            ind_tracked[i, 1] = np.argmin(np.abs(ind_tracked[i, 1] - YTemp))
            ind_tracked[i, 2] = np.argmin(np.abs(ind_tracked[i, 2] - ZTemp))

        outFile.write('    indTracked = reshape((/')
        outFile.write(','.join(map(str, ind_tracked.flatten())))
        outFile.write('/),shape(indTracked))\n')

    outFile.close()

def writeMatrix(outFile, name, var):
    outFile.write(f'    {name} = reshape((/')
    outFile.write(','.join(f'{v:.20f}d0' for v in var.flatten()))
    outFile.write(f'/),shape({name}))\n')



