import numpy as np

# Flow parameters
Re = 1000
Ma = 0.5
L = 4
H = 1

# Case name
caseName = f"Poiseulle-Re{Re}-Ma{str(Ma).replace('.', '')}-L{L}-H{H}"

# Domain decomposition
p_row = 10
p_col = 1

# Flow parameters
flowParameters = {
    'Re': Re,
    'Ma': Ma,
    'Pr': 0.71,
    'gamma': 1.4,
    'T0': 300
}

# Domain parameters
domain = {
    'xi': 0,
    'xf': L,
    'yi': -H/2,
    'yf': H/2,
    'zi': 0,
    'zf': 1
}

# Flow type
flowType = {
    'name': 'poiseulleFlow',
    'initial': {
        'type': 'poiseulle',  # uniform, blasius or file
        'addNoise': 1e-4
    }
}

# Mesh parameters
mesh = {
    'x': {
        'n': 100,
        'type': 'uniform',
        'matchFixed': 2,
        'periodic': True,
        'fixPeriodicDomainSize': False,
        'extraRefinement': 0,
        'buffer': {'i': {'n': 0}, 'f': {'n': 0}}
    },
    'y': {
        'n': 100,
        'type': 'tanh',
        'local': 'b',
        'par': 2,
        'matchFixed': 2,
        'periodic': False,
        'fixPeriodicDomainSize': False,
        'extraRefinement': 0,
        'buffer': {'i': {'n': 0}, 'f': {'n': 0}}
    },
    'z': {
        'n': 1,
        'type': 'uniform',
        'matchFixed': 2,
        'periodic': True,
        'fixPeriodicDomainSize': True,
        'extraRefinement': 0,
        'buffer': {'i': {'n': 0}, 'f': {'n': 0}}
    },
    'trackedPoints': [],
    'trackedNorm': True
}

# Time control
time = {
    'control': 'cfl',  # If 'dt', qtimes and tmax are in number of iterations and the step size is fixed to dt
                       # If 'cfl', qtimes and tmax are in non-dimensional time and dt defines the maximum step size
    'dt': 1,
    'maxCFL': 1.3,
    'qtimes': 1,
    'tmax': 1000
}

logAll = 25

# Numerical methods
numMethods = {
    'spatialDerivs': 'SL6',  # SL4 or EX2
    'spatialDerivsBuffer': 'EX4',
    'timeStepping': 'RK4',  # RK4 or Euler
    'neumannOrder': 6,
    'neumann2Order': 2,
    'spatialFilterStrength': 0.49,  # -0.5 < alpha < 0.5
    'spatialFilterTime': 0,  # Characteristic time of the spatial filter (optional, default = 0)
    'filterDirections': [1, 1, 1],
    'filterBorders': 'reducedOrder',
    'filterBordersStartX': False,
    'filterBordersEndX': False,
    'filterBordersEndY': False,
    'SFD': {
        'type': 0,  # 0 = off, 1 = whole domain, 2 = buffer zone only
        'X': 0.05,
        'Delta': 10,
        'applyY': False
    }
}

# Uncomment and modify if needed:
# numMethods['SFD']['extraRegion'] = [{
#     'location': [(x1+x2)/2, 0, 0],
#     'size': [L, D, float('inf')],
#     'X': 0.1
# }]


