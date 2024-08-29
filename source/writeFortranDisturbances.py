import os

# na versão original em matlab, disturbTypes é retornada como saída da função, mas ela não parece ser usada em qualquer outro lugar (gigiaero - 27/08/2024) 

def writeFortranDisturbances(case_name, bi, tridimensional):
    """
    This function creates the runDisturbances.F90 and disturbances.F90 files.
    runDisturbances.F90 is called when the boundary conditions are applied. It calls the appropriate routines for each domain slice.
    runForcings.F90 is called at the end of the Navier Stokes equations routine. It calls the appropriate routines for each domain slice.
    disturbances.F90 is a Fortran module file, which concatenates all the disturbances routines from the disturbances folder.
    """
    
    nProcs = len(bi)
    
    outFileDisturb = open(f"{case_name}/bin/runDisturbances.F90", 'w')
    if not tridimensional:
        outFileForcing = open(f"{case_name}/bin/runForcings2D.F90", 'w')
        os.system(f"touch {case_name}/bin/runForcings3D.F90")
    else:
        outFileForcing = open(f"{case_name}/bin/runForcings3D.F90", 'w')
        os.system(f"touch {case_name}/bin/runForcings2D.F90")
    
    outFileDisturb.write('    select case (nrank)\n')
    outFileForcing.write('    select case (nrank)\n')
    
    disturbTypes = []
    
    for i in range(nProcs):
        outFileDisturb.write(f'        case ({i})\n')
        outFileForcing.write(f'        case ({i})\n')
        
        for j in range(len(bi[i]['disturb'])):
            if bi[i]['disturb'][j]['forcing']:
                outFile = outFileForcing
            else:
                outFile = outFileDisturb
            
            di = bi[i]['disturb'][j]
            
            disturbTypes.append(di['type'])
            
            nx = di['ind'][1] - di['ind'][0] + 1
            ny = di['ind'][3] - di['ind'][2] + 1
            nz = di['ind'][5] - di['ind'][4] + 1
            
            outFile.write(f'            call {di["type"]}({nx},{ny},{nz},(/')
            
            for k in range(nx - 1):
                outFile.write(f'{di["X"][k]:.20f}d0,')
            outFile.write(f'{di["X"][-1]:.20f}d0/),(/')
            for k in range(ny - 1):
                outFile.write(f'{di["Y"][k]:.20f}d0,')
            outFile.write(f'{di["Y"][-1]:.20f}d0/),(/')
            for k in range(nz - 1):
                outFile.write(f'{di["Z"][k]:.20f}d0,')
            if bi[i]['disturb'][j]['forcing']:
                outFile.write(f'{di["Z"][-1]:.20f}d0/)')
            else:
                outFile.write(f'{di["Z"][-1]:.20f}d0/),t')
            
            if bi[i]['disturb'][j]['forcing']:
                for k in range(len(di['var'])):
                    outFile.write(f',{di["var"][k]}({di["ind"][0]}:{di["ind"][1]},{di["ind"][2]}:{di["ind"][3]},{di["ind"][4]}:{di["ind"][5]})')
                for k in range(len(di['var'])):
                    outFile.write(f',d{di["var"][k]}({di["ind"][0]}:{di["ind"][1]},{di["ind"][2]}:{di["ind"][3]},{di["ind"][4]}:{di["ind"][5]})')
            else:
                for k in range(len(di['var'])):
                    outFile.write(f',{di["var"][k]}({di["ind"][0]}:{di["ind"][1]},{di["ind"][2]}:{di["ind"][3]},{di["ind"][4]}:{di["ind"][5]})')
            
            for k in range(len(di['par'])):
                if isinstance(di['par'][k], (int, float)):
                    if isinstance(di['par'][k], (list, tuple)):
                        outFile.write(',(/')
                        for kk in range(len(di['par'][k]) - 1):
                            outFile.write(f'{di["par"][k][kk]:.20f}d0,')
                        outFile.write(f'{di["par"][k][-1]:.20f}d0/)')
                    else:
                        outFile.write(f',{di["par"][k]:.20f}d0')
                else:
                    outFile.write(f",'{di['par'][k]}'")
            outFile.write(')\n')
    
    outFileForcing.write('    end select\n')
    outFileDisturb.write('    end select\n')
    
    outFileForcing.close()
    outFileDisturb.close()
    
    disturbTypes = list(set(disturbTypes))
    
    outFile = open(f"{case_name}/bin/disturbances.F90", 'w')
    outFile.write('    module disturbances\n\n    contains\n\n')
    
    for disturb_type in disturbTypes:
        if disturb_type != 'holdInlet':
            source_file = open(f'source/disturbances/{disturb_type}.F90', 'r')
        else:
            source_file = open(f'source/Fortran/{disturb_type}.F90', 'r')
        line = source_file.readline()
        while line:
            outFile.write(f'{line}\n')
            line = source_file.readline()
        source_file.close()
        
        outFile.write('\n')
    
    outFile.write('\n    end module\n')
    outFile.close()



