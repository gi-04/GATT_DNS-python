import numpy as np
import os

def get_boundary_conditions(flowType, mesh, flowParameters, neumannOrder):
    # This function runs the relevant routine for the flowType defined in the parameters and outputs the variables in a structure
    # The disturbances are also added to the structure
    # The stencil and coefficients for Neumann conditions is defined here

    ## Define base value for gamma, P and E
    gamma = flowParameters['gamma']
    E0 = 1 / ((gamma**2 - gamma) * flowParameters['Ma']**2)
    P0 = (gamma - 1) * E0

    ## Create variables that will receive boundary information
    var = []  # Variable to which the condition is applied
    type = []  # Type of boundary. Either 'dir' for a Dirichlet condition or 'neu' for a Neumann condition
    dir = []  # Direction of the derivative for the Neumann condition (xi, xf, yi, yf, zi or zf)
    val = []  # Value the Variable assumes in a Dirichlet condition or value of the derivative in a Neumann condition
    xi = []  # Starting j index for the condition
    xf = []  # Ending j index for the variable
    yi = []  # Starting i index for the condition
    yf = []  # Ending i index for the condition
    zi = []  # Starting k index for the condition
    zf = []  # Ending k index for the condition

    ## Check the type of boundaries and call the appropriate subroutine
    boundary_file = os.path.join('source', 'boundaries', f"{flowType['name']}.py")
    if os.path.exists(boundary_file):
        exec(open(boundary_file).read())
    else:
        raise FileNotFoundError('Boundary condition file not found, check source/boundaries folder')

    ## Join all relevant variables in the output structure
    boundary = {
        'var': var,
        'type': type,
        'dir': dir,
        'val': val,
        'xi': xi,
        'xf': xf,
        'yi': yi,
        'yf': yf,
        'zi': zi,
        'zf': zf,
        'flowRegion': flowRegion,
        'E0': E0,
        'wall': {
            'up': wallUpLimits,
            'down': wallDownLimits,
            'front': wallFrontLimits,
            'back': wallBackLimits,
            'right': wallRightLimits,
            'left': wallLeftLimits
        },
        'corners': corners,
        'gamma': gamma
    }

    ## Get wall region
    boundary['insideWall'] = ~flowRegion
    for i in range(1, 7):
        currentWall = {
            1: wallUpLimits,
            2: wallDownLimits,
            3: wallFrontLimits,
            4: wallBackLimits,
            5: wallRightLimits,
            6: wallLeftLimits
        }[i]
        for j in range(currentWall.shape[0]):
            boundary['insideWall'][currentWall[j, 0]:currentWall[j, 1], currentWall[j, 2]:currentWall[j, 3], currentWall[j, 4]:currentWall[j, 5]] = False

    ## Get coefficients for Neumann boundary conditions
    NC1 = [
        [-1, 1], [-3/2, 2, -1/2], [-11/6, 3, -3/2, 1/3], 
        [-25/12, 4, -3, 4/3, -1/4], [-137/60, 5, -5, 10/3, -5/4, 1/5], 
        [-49/20, 6, -15/2, 20/3, -15/4, 6/5, -1/6]
    ]

    NC2 = [
        [1, -2, 1], [2, -5, 4, -1], [35/12, -26/3, 19/2, -14/3, 11/12], 
        [15/4, -77/6, 107/6, -13, 61/12, -5/6], [203/45, -87/5, 117/4, -254/9, 33/2, -27/5, 137/180], 
        [469/90, -223/10, 879/20, -949/18, 41, -201/10, 1019/180, -7/10]
    ]

    boundary['neumannCoeffs'] = -np.array(NC1[neumannOrder[0]][1:]) / NC1[neumannOrder[0]][0]
    boundary['neumann2Coeffs'] = -np.array(NC2[neumannOrder[1]][1:]) / NC2[neumannOrder[1]][0]

    ## Add all disturbances to boundary structure
    boundary['disturb'] = []
    if 'disturb' in flowType:
        for i, disturb in enumerate(flowType['disturb']):
            if disturb['active']:
                disturb_info = {
                    'type': disturb['type'],
                    'forcing': disturb.get('forcing', False),
                    'par': disturb.get('par', [])
                }

                if isinstance(disturb_info['par'], np.ndarray):
                    disturb_info['par'] = disturb_info['par'].tolist()

                disturb_info['var'] = disturb['var']

                xi, xf = disturb['x']
                xi = mesh['X'][0] if np.isinf(xi) and xi < 0 else xi
                xf = mesh['X'][0] if np.isinf(xf) and xf < 0 else xf
                xi = mesh['X'][-1] if np.isinf(xi) and xi > 0 else xi
                xf = mesh['X'][-1] if np.isinf(xf) and xf > 0 else xf

                yi, yf = disturb['y']
                yi = mesh['Y'][0] if np.isinf(yi) and yi < 0 else yi
                yf = mesh['Y'][0] if np.isinf(yf) and yf < 0 else yf
                yi = mesh['Y'][-1] if np.isinf(yi) and yi > 0 else yi
                yf = mesh['Y'][-1] if np.isinf(yf) and yf > 0 else yf

                zi, zf = disturb['z']
                zi = mesh['Z'][0] if np.isinf(zi) and zi < 0 else zi
                zf = mesh['Z'][0] if np.isinf(zf) and zf < 0 else zf
                zi = mesh['Z'][-1] if np.isinf(zi) and zi > 0 else zi
                zf = mesh['Z'][-1] if np.isinf(zf) and zf > 0 else zf

                disturb_info['extraNodes'] = disturb.get('extraNodes', [0, 0, 0, 0, 0, 0])

                disturb_info['ind'] = [
                    np.searchsorted(mesh['X'], xi, side='left') - disturb_info['extraNodes'][0],
                    np.searchsorted(mesh['X'], xf, side='right') + disturb_info['extraNodes'][1],
                    np.searchsorted(mesh['Y'], yi, side='left') - disturb_info['extraNodes'][2],
                    np.searchsorted(mesh['Y'], yf, side='right') + disturb_info['extraNodes'][3],
                    np.searchsorted(mesh['Z'], zi, side='left') - disturb_info['extraNodes'][4],
                    np.searchsorted(mesh['Z'], zf, side='right') + disturb_info['extraNodes'][5]
                ]

                boundary['disturb'].append(disturb_info)

    return boundary, mesh



