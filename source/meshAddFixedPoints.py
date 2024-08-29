# This script creates a list of interest points in the mesh. If mesh.matchFixed is true, the mesh will be slightly transformed to include these points

import numpy as np
# from __main__ import *
# from parameters import *

# transformei este script em uma função (gigiaero - 27/08/2024)

def meshAddFixedPoints(mesh,flowType,domain):

    # mesh = {'x': {},
    #         'y': {},
    #         'z': {}}
    mesh['x']['fixPoints'] = [0]
    mesh['y']['fixPoints'] = [0]
    mesh['z']['fixPoints'] = []

    if 'cav' in flowType:
        for cav in flowType['cav']:
            mesh['x']['fixPoints'].extend(cav['x'])
            mesh['y']['fixPoints'].extend(cav['y'])
            mesh['z']['fixPoints'].extend(cav['z'])

    if 'rug' in flowType:
        for rug in flowType['rug']:
            mesh['x']['fixPoints'].extend(rug['x'])
            mesh['y']['fixPoints'].extend(rug['y'])
            mesh['z']['fixPoints'].extend(rug['z'])

    if 'disturb' in flowType:
        for disturb in flowType['disturb']:
            if 'fitPoints' in disturb and disturb['fitPoints']:
                mesh['x']['fixPoints'].extend(disturb['x'])
                mesh['y']['fixPoints'].extend(disturb['y'])
                mesh['z']['fixPoints'].extend(disturb['z'])
            else:
                if 'fitPointsX' in disturb and disturb['fitPointsX']:
                    mesh['x']['fixPoints'].extend(disturb['x'])
                if 'fitPointsY' in disturb and disturb['fitPointsY']:
                    mesh['y']['fixPoints'].extend(disturb['y'])
                if 'fitPointsZ' in disturb and disturb['fitPointsZ']:
                    mesh['z']['fixPoints'].extend(disturb['z'])

    if 'trackedPoints' in mesh and 'fitTrackedPoints' in mesh and mesh['fitTrackedPoints']:
        for point in mesh['trackedPoints']:
            mesh['x']['fixPoints'].extend(point[0])
            mesh['y']['fixPoints'].extend(point[1])
            mesh['z']['fixPoints'].extend(point[2])

    # Remove duplicates and sort the fixPoints
    # mesh['x']['fixPoints'] = sorted(set(mesh['x']['fixPoints'])) # conversões originais que não funcionaram - ainda falta testar (gigiaero - 26/08/2024)
    # mesh['y']['fixPoints'] = sorted(set(mesh['y']['fixPoints']))
    # mesh['z']['fixPoints'] = sorted(set(mesh['z']['fixPoints']))
    mesh['x']['fixPoints'] = np.unique(mesh['x']['fixPoints'])
    mesh['y']['fixPoints'] = np.unique(mesh['y']['fixPoints'])
    mesh['z']['fixPoints'] = np.unique(mesh['z']['fixPoints'])

    # Remove points that are out of the domain or are infinite
    mesh['x']['fixPoints'] = [x for x in mesh['x']['fixPoints'] if not np.isinf(x) and domain['xi'] <= x <= domain['xf']]
    mesh['y']['fixPoints'] = [y for y in mesh['y']['fixPoints'] if not np.isinf(y) and domain['yi'] <= y <= domain['yf']]
    mesh['z']['fixPoints'] = [z for z in mesh['z']['fixPoints'] if not np.isinf(z) and domain['zi'] <= z <= domain['zf']]

    return mesh
