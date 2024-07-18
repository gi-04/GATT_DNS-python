import os

def write_fortran_disturbances(case_name, bi, tridimensional):
    """
    This function creates the runDisturbances.F90 and disturbances.F90 files.
    runDisturbances.F90 is called when the boundary conditions are applied. It calls the appropriate routines for each domain slice.
    runForcings.F90 is called at the end of the Navier Stokes equations routine. It calls the appropriate routines for each domain slice.
    disturbances.F90 is a Fortran module file, which concatenates all the disturbances routines from the disturbances folder.
    """
    
    n_procs = len(bi)
    
    out_file_disturb = open(f"{case_name}/bin/runDisturbances.F90", 'w')
    if not tridimensional:
        out_file_forcing = open(f"{case_name}/bin/runForcings2D.F90", 'w')
        os.system(f"touch {case_name}/bin/runForcings3D.F90")
    else:
        out_file_forcing = open(f"{case_name}/bin/runForcings3D.F90", 'w')
        os.system(f"touch {case_name}/bin/runForcings2D.F90")
    
    out_file_disturb.write('    select case (nrank)\n')
    out_file_forcing.write('    select case (nrank)\n')
    
    disturb_types = []
    
    for i in range(n_procs):
        out_file_disturb.write(f'        case ({i})\n')
        out_file_forcing.write(f'        case ({i})\n')
        
        for j in range(len(bi[i]['disturb'])):
            if bi[i]['disturb'][j]['forcing']:
                out_file = out_file_forcing
            else:
                out_file = out_file_disturb
            
            di = bi[i]['disturb'][j]
            
            disturb_types.append(di['type'])
            
            nx = di['ind'][1] - di['ind'][0] + 1
            ny = di['ind'][3] - di['ind'][2] + 1
            nz = di['ind'][5] - di['ind'][4] + 1
            
            out_file.write(f'            call {di["type"]}({nx},{ny},{nz},(/')
            
            for k in range(nx - 1):
                out_file.write(f'{di["X"][k]:.20f}d0,')
            out_file.write(f'{di["X"][-1]:.20f}d0/),(/')
            for k in range(ny - 1):
                out_file.write(f'{di["Y"][k]:.20f}d0,')
            out_file.write(f'{di["Y"][-1]:.20f}d0/),(/')
            for k in range(nz - 1):
                out_file.write(f'{di["Z"][k]:.20f}d0,')
            if bi[i]['disturb'][j]['forcing']:
                out_file.write(f'{di["Z"][-1]:.20f}d0/)')
            else:
                out_file.write(f'{di["Z"][-1]:.20f}d0/),t')
            
            if bi[i]['disturb'][j]['forcing']:
                for k in range(len(di['var'])):
                    out_file.write(f',{di["var"][k]}({di["ind"][0]}:{di["ind"][1]},{di["ind"][2]}:{di["ind"][3]},{di["ind"][4]}:{di["ind"][5]})')
                for k in range(len(di['var'])):
                    out_file.write(f',d{di["var"][k]}({di["ind"][0]}:{di["ind"][1]},{di["ind"][2]}:{di["ind"][3]},{di["ind"][4]}:{di["ind"][5]})')
            else:
                for k in range(len(di['var'])):
                    out_file.write(f',{di["var"][k]}({di["ind"][0]}:{di["ind"][1]},{di["ind"][2]}:{di["ind"][3]},{di["ind"][4]}:{di["ind"][5]})')
            
            for k in range(len(di['par'])):
                if isinstance(di['par'][k], (int, float)):
                    if isinstance(di['par'][k], (list, tuple)):
                        out_file.write(',(/')
                        for kk in range(len(di['par'][k]) - 1):
                            out_file.write(f'{di["par"][k][kk]:.20f}d0,')
                        out_file.write(f'{di["par"][k][-1]:.20f}d0/)')
                    else:
                        out_file.write(f',{di["par"][k]:.20f}d0')
                else:
                    out_file.write(f",'{di['par'][k]}'")
            out_file.write(')\n')
    
    out_file_forcing.write('    end select\n')
    out_file_disturb.write('    end select\n')
    
    out_file_forcing.close()
    out_file_disturb.close()
    
    disturb_types = list(set(disturb_types))
    
    out_file = open(f"{case_name}/bin/disturbances.F90", 'w')
    out_file.write('    module disturbances\n\n    contains\n\n')
    
    for disturb_type in disturb_types:
        if disturb_type != 'holdInlet':
            source_file = open(f'source/disturbances/{disturb_type}.F90', 'r')
        else:
            source_file = open(f'source/Fortran/{disturb_type}.F90', 'r')
        line = source_file.readline()
        while line:
            out_file.write(f'{line}\n')
            line = source_file.readline()
        source_file.close()
        
        out_file.write('\n')
    
    out_file.write('\n    end module\n')
    out_file.close()



