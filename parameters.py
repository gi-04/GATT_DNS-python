import numpy as np

delta = 1
Red = 600
Ma = 0.5
L = 10
D = 5
x1 = delta**2 * Red / (1.72**2 * delta**2)
x2 = x1 + L
xEnd = x2 + 200
delta99end = 5 * xEnd / np.sqrt(Red * xEnd)
yEnd = 5 * delta99end

# Case name
caseName = f'Red{Red}-Ma{str(Ma).replace(".", "")}-L{L}-D{D}'

# Domain decomposition
p_row = 10
p_col = 1

# Flow parameters
flowParameters = {
    'Re': Red,
    'Ma': Ma,
    'Pr': 0.71,
    'gamma': 1.4,
    'T0': 300
}

# Domain parameters
domain = {
    'xi': 0,
    'xf': xEnd,
    'yi': -D,
    'yf': yEnd,
    'zi': 0,
    'zf': 1
}

# Flow type
flowType = {
    'name': 'boundaryLayerIsothermal',
    'initial': {
        'type': 'blasius'
    },
    'cav': [{
        'x': [x1, x2],
        'y': [-D, 0],
        'z': [-np.inf, np.inf]
    }],
    'disturb': [{
        'x': [25, 50],
        'y': [0, 0],
        'z': [-np.inf, np.inf],
        'var': 'V',
        'type': 'packet_2d',
        'extraNodes': [0, 0, 0, 0, 0, 0],
        'par': [0.02, 50, 1e-5],
        'active': False,
        'fitPoints': False
    }]
}

# Mesh parameters
mesh = {
    'x': {
        'd0': 4,
        'type': 'attractors',
        'attractorPoints': [],
        'attractorStrength': [],
        'attractorSize': [],
        'attractorRegions': [
            [x1-100, x1-50, x2+50, x2+100, 3],
            [x1-20, x1, x2+10, x2+20, 16],
            [x2-5, x2-0.1, x2, x2+2, 30],
            [x1-2, x1, x1+0.1, x1+5, 20],
            [-30, -10, 10, 30, 1]
        ],
        'matchFixed': 2,
        'periodic': False,
        'fixPeriodicDomainSize': False,
        'extraRefinement': 0,
        'buffer': {
            'i': {
                'n': 40,
                'type': 'sigmoid',
                'stretching': 30,
                'transition': 0.2
            },
            'f': {
                'n': 30,
                'type': 'sigmoid',
                'stretching': 15,
                'transition': 0.2
            }
        }
    },
    'y': {
        'd0': 1,
        'type': 'attractors',
        'attractorPoints': [],
        'attractorStrength': [],
        'attractorSize': [],
        'attractorRegions': [[-2*D, -D/10, D/10, 2*delta99end, 12]],
        'matchFixed': 2,
        'periodic': False,
        'fixPeriodicDomainSize': False,
        'extraRefinement': 0,
        'buffer': {
            'i': {
                'n': 0
            },
            'f': {
                'n': 40,
                'type': 'sigmoid',
                'stretching': 30,
                'transition': 0.2
            }
        }
    },
    'z': {
        'n': 1,
        'type': 'uniform',
        'matchFixed': 2,
        'periodic': True,
        'fixPeriodicDomainSize': True,
        'extraRefinement': 0,
        'buffer': {
            'i': {
                'n': 0
            },
            'f': {
                'n': 0
            }
        }
    }
}

# Tracked points
trackedX = np.arange(0, domain['xf'] + 1, 50)
nProbes = len(trackedX)
mesh['trackedPoints'] = np.vstack((
    [37.5, 0, 0],
    [(x1 + 3 * x2) / 4, 0, 0],
    np.column_stack((trackedX, np.ones(nProbes), np.zeros(nProbes)))
))
mesh['trackedNorm'] = True

# Time control
time = {
    'control': 'cfl',
    'dt': 1,
    'maxCFL': 1.3,
    'qtimes': 100,
    'tmax': 10000
}

logAll = 25

# Numerical methods
numMethods = {
    'spatialDerivs': 'SL6',
    'spatialDerivsBuffer': 'EX4',
    'timeStepping': 'RK4',
    'neumannOrder': 6,
    'neumann2Order': 2,
    'spatialFilterStrength': 0.49,
    'spatialFilterTime': 0,
    'filterDirections': [1, 1, 1],
    'filterBorders': 'reducedOrder',
    'filterBordersStartX': False,
    'filterBordersEndX': False,
    'filterBordersEndY': False,
    'SFD': {
        'type': 2,
        'X': 0.05,
        'Delta': 10,
        'applyY': False
    }
}



