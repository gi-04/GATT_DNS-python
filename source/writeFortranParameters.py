import numpy as np

def writeFortranParameters(caseName, mesh, flowParameters, time, numMethods, logAll, p_row, p_col):
    """
    This function writes the parameters.F90 file, which contains basic parameters that will be used in the simulation.
    """
    with open(f'{caseName}/bin/parameters.F90', 'w') as outFile:
        outFile.write(f'    integer :: nx = {mesh["nx"]}\n')
        outFile.write(f'    integer :: ny = {mesh["ny"]}\n')
        outFile.write(f'    integer :: nz = {mesh["nz"]}\n\n')
        
        outFile.write(f'    real*8 :: Re = {flowParameters["Re"]:.20f}d0\n')
        outFile.write(f'    real*8 :: Ma = {flowParameters["Ma"]:.20f}d0\n')
        outFile.write(f'    real*8 :: Pr = {flowParameters["Pr"]:.20f}d0\n')
        outFile.write(f'    real*8 :: T0 = {flowParameters["T0"]:.20f}d0\n')
        outFile.write(f'    real*8 :: gamma = {flowParameters["gamma"]:.20f}d0\n\n')
        
        outFile.write(f'    real*8 :: dtmax = {time["dt"]:.20f}d0\n')
        outFile.write(f'    real*8 :: maxCFL = {time["maxCFL"]:.20f}d0\n')
        
        if logAll == 0:
            logAll = 2147483647
        outFile.write(f'    integer :: logAll = {logAll}\n')
        
        outFile.write(f'    integer :: nSave = {time["nStep"]}\n')
        
        if 'trackedNorm' not in mesh or not mesh['trackedNorm']:
            outFile.write('    real*8 :: trackedNorm = 0.d0\n')
        else:
            tracked_norm = 1 / ((flowParameters["gamma"]**2 - flowParameters["gamma"]) * flowParameters["Ma"]**2)
            outFile.write(f'    real*8 :: trackedNorm = {tracked_norm:.20f}d0\n')
        
        if time["control"] == 'dt':
            outFile.write('    integer :: timeControl = 1\n')
            outFile.write(f'    integer :: qTimesInt = {time["qtimes"]}\n')
            outFile.write('    real*8  :: qTimesReal\n\n')
            outFile.write(f'    integer :: tmaxInt = {time["tmax"]}\n')
            outFile.write('    real*8  :: tmaxReal\n\n')
        elif time["control"] == 'cfl':
            outFile.write('    integer :: timeControl = 2\n')
            outFile.write('    integer :: qTimesInt\n')
            outFile.write(f'    real*8  :: qtimesReal = {time["qtimes"]:.20f}d0\n\n')
            outFile.write('    integer :: tmaxInt\n')
            outFile.write(f'    real*8  :: tmaxReal = {time["tmax"]:.20f}d0\n\n')
        else:
            raise ValueError('Unrecognized type of time control. Use either dt or cfl')
        
        if numMethods["timeStepping"] == 'RK4':
            outFile.write('    integer :: timeStepping = 1\n')
        elif numMethods["timeStepping"] == 'Euler':
            outFile.write('    integer :: timeStepping = 2\n')
        elif numMethods["timeStepping"] == 'SSPRK3':
            outFile.write('    integer :: timeStepping = 3\n')
        else:
            raise ValueError('Unrecognized time stepping method')
        
        if 'SFD' in numMethods:
            outFile.write(f'    integer :: SFD = {numMethods["SFD"]["type"]}\n')
            outFile.write(f'    real*8 :: SFD_Delta = {numMethods["SFD"]["Delta"]:.20f}d0\n')
            outFile.write(f'    real*8 :: SFD_X_val = {numMethods["SFD"]["X"]:.20f}d0\n')
            outFile.write(f'    integer :: resumeMeanFlow = {numMethods["SFD"]["resume"]}\n\n')
        else:
            outFile.write('    integer :: SFD = 0\n')
            outFile.write('    real*8 :: SFD_Delta = 0.00000000000000000000d0\n')
            outFile.write('    real*8 :: SFD_X_val = 0.00000000000000000000d0\n')
            outFile.write('    integer :: resumeMeanFlow = 0\n\n')
        
        if 'spatialFilterTime' not in numMethods or numMethods['spatialFilterTime'] <= 0:
            outFile.write('    real*8 :: FilterCharTime = -1.00000000000000000000d0\n')
        else:
            outFile.write(f'    real*8 :: FilterCharTime = {numMethods["spatialFilterTime"]:.20f}d0\n')
        
        if mesh["nz"] == 1:
            dxmin = [1 / np.min(np.diff(mesh["X"])), 1 / np.min(np.diff(mesh["Y"])), 0]
        else:
            dxmin = [1 / np.min(np.diff(mesh["X"])), 1 / np.min(np.diff(mesh["Y"])), 1 / np.min(np.diff(mesh["Z"]))]
        
        if 'CFLignoreZ' in time and time['CFLignoreZ']:
            dxmin[2] = 0
        
        outFile.write(f'    real*8,dimension(3) :: dxmin = (/{dxmin[0]:.20f}d0,{dxmin[1]:.20f}d0,{dxmin[2]:.20f}d0/)\n\n')
        
        outFile.write(f'    integer :: p_row = {p_row}\n')
        outFile.write(f'    integer :: p_col = {p_col}\n')



