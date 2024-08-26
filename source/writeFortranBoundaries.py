def writeFortranBoundaries(case_name, bi):
    """
    This function writes the boundaryInfo.F90 file, it has the data for all boundary conditions for each domain slice
    """
    import numpy as np

    nProcs = len(bi)

    with open(f'{case_name}/bin/boundaryInfo.F90', 'w') as outFile:
        outFile.write('    select case (nrank)\n')

        vars = ['U', 'V', 'W', 'E', 'P']

        for i in range(nProcs):
            outFile.write(f'        case ({i})\n')

            # Dirichlet
            for j, var in enumerate(vars):
                if j == 0:
                    n = bi[i]['nUd']
                    ind = bi[i]['iUd']
                    val = bi[i]['vUd']
                elif j == 1:
                    n = bi[i]['nVd']
                    ind = bi[i]['iVd']
                    val = bi[i]['vVd']
                elif j == 2:
                    n = bi[i]['nWd']
                    ind = bi[i]['iWd']
                    val = bi[i]['vWd']
                elif j == 3:
                    n = bi[i]['nEd']
                    ind = bi[i]['iEd']
                    val = bi[i]['vEd']
                elif j == 4:
                    n = bi[i]['nPd']
                    ind = bi[i]['iPd']
                    val = bi[i]['vPd']

                outFile.write(f'            n{var}d = {n}\n')

                if n > 0:
                    outFile.write(f'            allocate(i{var}d({n},6))\n')
                    outFile.write(f'            allocate(v{var}d({n}))\n')

                    outFile.write(f'            i{var}d = reshape((/')
                    outFile.write(','.join(map(str, ind[:-1])) + f',{ind[-1]}/),shape(i{var}d))\n')

                    outFile.write(f'            v{var}d = (/')
                    outFile.write(','.join(f'{v:.20f}d0' for v in val[:-1]) + f',{val[-1]:.20f}d0/)\n\n')

            # Neumann
            for j, var in enumerate(vars):
                if j == 0:
                    n = bi[i]['nUn']
                    ind = bi[i]['iUn']
                    dir = bi[i]['dUn']
                elif j == 1:
                    n = bi[i]['nVn']
                    ind = bi[i]['iVn']
                    dir = bi[i]['dVn']
                elif j == 2:
                    n = bi[i]['nWn']
                    ind = bi[i]['iWn']
                    dir = bi[i]['dWn']
                elif j == 3:
                    n = bi[i]['nEn']
                    ind = bi[i]['iEn']
                    dir = bi[i]['dEn']
                elif j == 4:
                    n = bi[i]['nPn']
                    ind = bi[i]['iPn']
                    dir = bi[i]['dPn']

                outFile.write(f'            n{var}n = {n}\n')

                if n > 0:
                    outFile.write(f'            allocate(i{var}n({n},6))\n')
                    outFile.write(f'            allocate(d{var}n({n}))\n')

                    outFile.write(f'            i{var}n = reshape((/')
                    outFile.write(','.join(map(str, ind[:-1])) + f',{ind[-1]}/),shape(i{var}n))\n')

                    outFile.write(f'            d{var}n = (/')
                    outFile.write(','.join(map(str, dir[:-1])) + f',{dir[-1]}/)\n\n')
                else:
                    outFile.write(f'            allocate(i{var}n(1,6))\n')
                    outFile.write(f'            allocate(d{var}n(1))\n\n')

            # Second derivative
            for j, var in enumerate(vars):
                if j == 0:
                    n = bi[i]['nUs']
                    ind = bi[i]['iUs']
                    dir = bi[i]['dUs']
                elif j == 1:
                    n = bi[i]['nVs']
                    ind = bi[i]['iVs']
                    dir = bi[i]['dVs']
                elif j == 2:
                    n = bi[i]['nWs']
                    ind = bi[i]['iWs']
                    dir = bi[i]['dWs']
                elif j == 3:
                    n = bi[i]['nEs']
                    ind = bi[i]['iEs']
                    dir = bi[i]['dEs']
                elif j == 4:
                    n = bi[i]['nPs']
                    ind = bi[i]['iPs']
                    dir = bi[i]['dPs']

                outFile.write(f'            n{var}s = {n}\n')

                if n > 0:
                    outFile.write(f'            allocate(i{var}s({n},6))\n')
                    outFile.write(f'            allocate(d{var}s({n}))\n')

                    outFile.write(f'            i{var}s = reshape((/')
                    outFile.write(','.join(map(str, ind[:-1])) + f',{ind[-1]}/),shape(i{var}s))\n')

                    outFile.write(f'            d{var}s = (/')
                    outFile.write(','.join(map(str, dir[:-1])) + f',{dir[-1]}/)\n\n')
                else:
                    outFile.write(f'            allocate(i{var}s(1,6))\n')
                    outFile.write(f'            allocate(d{var}s(1))\n\n')

            # Corners
            outFile.write(f'            cN = {bi[i]["cN"]}\n')
            outFile.write(f'            allocate(cL({max(bi[i]["cN"], 1)},6))\n')
            outFile.write(f'            allocate(cD({max(bi[i]["cN"], 1)},3))\n')
            outFile.write(f'            allocate(cAdiabatic({max(bi[i]["cN"], 1)}))\n')

            if bi[i]['cN'] > 0:
                outFile.write(f'            cL = reshape((/')
                outFile.write(','.join(map(str, bi[i]['cL'][:-1])) + f',{bi[i]["cL"][-1]}/),shape(cL))\n')

                outFile.write(f'            cD = reshape((/')
                outFile.write(','.join(map(str, bi[i]['cD'][:-1])) + f',{bi[i]["cD"][-1]}/),shape(cD))\n')

                outFile.write(f'            cAdiabatic = (/')
                outFile.write(','.join(map(str, bi[i]['adiabatic'][:-1])) + f',{bi[i]["adiabatic"][-1]}/)\n\n')

        outFile.write('    end select\n')



