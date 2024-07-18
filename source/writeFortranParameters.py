import numpy as np

def write_fortran_parameters(case_name, mesh, flow_parameters, time, num_methods, log_all, p_row, p_col):
    """
    This function writes the parameters.F90 file, which contains basic parameters that will be used in the simulation.
    """
    with open(f'{case_name}/bin/parameters.F90', 'w') as out_file:
        out_file.write(f'    integer :: nx = {mesh["nx"]}\n')
        out_file.write(f'    integer :: ny = {mesh["ny"]}\n')
        out_file.write(f'    integer :: nz = {mesh["nz"]}\n\n')
        
        out_file.write(f'    real*8 :: Re = {flow_parameters["Re"]:.20f}d0\n')
        out_file.write(f'    real*8 :: Ma = {flow_parameters["Ma"]:.20f}d0\n')
        out_file.write(f'    real*8 :: Pr = {flow_parameters["Pr"]:.20f}d0\n')
        out_file.write(f'    real*8 :: T0 = {flow_parameters["T0"]:.20f}d0\n')
        out_file.write(f'    real*8 :: gamma = {flow_parameters["gamma"]:.20f}d0\n\n')
        
        out_file.write(f'    real*8 :: dtmax = {time["dt"]:.20f}d0\n')
        out_file.write(f'    real*8 :: maxCFL = {time["maxCFL"]:.20f}d0\n')
        
        if log_all == 0:
            log_all = 2147483647
        out_file.write(f'    integer :: logAll = {log_all}\n')
        
        out_file.write(f'    integer :: nSave = {time["nStep"]}\n')
        
        if 'trackedNorm' not in mesh or not mesh['trackedNorm']:
            out_file.write('    real*8 :: trackedNorm = 0.d0\n')
        else:
            tracked_norm = 1 / ((flow_parameters["gamma"]**2 - flow_parameters["gamma"]) * flow_parameters["Ma"]**2)
            out_file.write(f'    real*8 :: trackedNorm = {tracked_norm:.20f}d0\n')
        
        if time["control"] == 'dt':
            out_file.write('    integer :: timeControl = 1\n')
            out_file.write(f'    integer :: qTimesInt = {time["qtimes"]}\n')
            out_file.write('    real*8  :: qTimesReal\n\n')
            out_file.write(f'    integer :: tmaxInt = {time["tmax"]}\n')
            out_file.write('    real*8  :: tmaxReal\n\n')
        elif time["control"] == 'cfl':
            out_file.write('    integer :: timeControl = 2\n')
            out_file.write('    integer :: qTimesInt\n')
            out_file.write(f'    real*8  :: qtimesReal = {time["qtimes"]:.20f}d0\n\n')
            out_file.write('    integer :: tmaxInt\n')
            out_file.write(f'    real*8  :: tmaxReal = {time["tmax"]:.20f}d0\n\n')
        else:
            raise ValueError('Unrecognized type of time control. Use either dt or cfl')
        
        if num_methods["timeStepping"] == 'RK4':
            out_file.write('    integer :: timeStepping = 1\n')
        elif num_methods["timeStepping"] == 'Euler':
            out_file.write('    integer :: timeStepping = 2\n')
        elif num_methods["timeStepping"] == 'SSPRK3':
            out_file.write('    integer :: timeStepping = 3\n')
        else:
            raise ValueError('Unrecognized time stepping method')
        
        if 'SFD' in num_methods:
            out_file.write(f'    integer :: SFD = {num_methods["SFD"]["type"]}\n')
            out_file.write(f'    real*8 :: SFD_Delta = {num_methods["SFD"]["Delta"]:.20f}d0\n')
            out_file.write(f'    real*8 :: SFD_X_val = {num_methods["SFD"]["X"]:.20f}d0\n')
            out_file.write(f'    integer :: resumeMeanFlow = {num_methods["SFD"]["resume"]}\n\n')
        else:
            out_file.write('    integer :: SFD = 0\n')
            out_file.write('    real*8 :: SFD_Delta = 0.00000000000000000000d0\n')
            out_file.write('    real*8 :: SFD_X_val = 0.00000000000000000000d0\n')
            out_file.write('    integer :: resumeMeanFlow = 0\n\n')
        
        if 'spatialFilterTime' not in num_methods or num_methods['spatialFilterTime'] <= 0:
            out_file.write('    real*8 :: FilterCharTime = -1.00000000000000000000d0\n')
        else:
            out_file.write(f'    real*8 :: FilterCharTime = {num_methods["spatialFilterTime"]:.20f}d0\n')
        
        if mesh["nz"] == 1:
            dxmin = [1 / np.min(np.diff(mesh["X"])), 1 / np.min(np.diff(mesh["Y"])), 0]
        else:
            dxmin = [1 / np.min(np.diff(mesh["X"])), 1 / np.min(np.diff(mesh["Y"])), 1 / np.min(np.diff(mesh["Z"]))]
        
        if 'CFLignoreZ' in time and time['CFLignoreZ']:
            dxmin[2] = 0
        
        out_file.write(f'    real*8,dimension(3) :: dxmin = (/{dxmin[0]:.20f}d0,{dxmin[1]:.20f}d0,{dxmin[2]:.20f}d0/)\n\n')
        
        out_file.write(f'    integer :: p_row = {p_row}\n')
        out_file.write(f'    integer :: p_col = {p_col}\n')



