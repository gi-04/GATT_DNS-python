import numpy as np

# Flow parameters
Re = 1000
Ma = 0.5
L = 1
D = 1

# Case name
case_name = f"LidDriven-Re{Re}-Ma{str(Ma).replace('.', '')}-L{L}-D{D}"

# Domain decomposition
p_row = 10
p_col = 1

# Flow parameters
flow_parameters = {
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
    'yi': 0,
    'yf': D,
    'zi': 0,
    'zf': 1
}

# Flow type
flow_type = {
    'name': 'lidDrivenFlow',
    'initial': {
        'type': 'uniform',  # uniform, blasius or file
        'U0': 0
    }
}

# Mesh parameters
mesh = {
    'x': {'n': 100, 'type': 'tanh', 'local': 'b', 'par': 2, 'matchFixed': 2, 'periodic': False, 'fixPeriodicDomainSize': False, 'extraRefinement': 0},
    'y': {'n': 100, 'type': 'tanh', 'local': 'b', 'par': 2, 'matchFixed': 2, 'periodic': False, 'fixPeriodicDomainSize': False, 'extraRefinement': 0},
    'z': {'n': 1, 'type': 'uniform', 'matchFixed': 2, 'periodic': True, 'fixPeriodicDomainSize': True, 'extraRefinement': 0}
}

for axis in ['x', 'y', 'z']:
    mesh[axis]['buffer'] = {'i': {'n': 0}, 'f': {'n': 0}}

# Add lid driven movement
flow_type['disturb'] = [{
    'x': [0, L],
    'y': [D, D],
    'z': [-np.inf, np.inf],
    'var': 'U',
    'type': 'lidDrivenMovement',
    'extraNodes': [0, 0, 0, 0, 0, 0],
    'par': [1],  # U0
    'active': True,
    'fitPoints': True
}]

# Time control
time = {
    'control': 'cfl',
    'dt': 1,
    'maxCFL': 1.3,
    'qtimes': 1,
    'tmax': 1000
}

log_all = 25

mesh['trackedPoints'] = []
mesh['trackedNorm'] = True

# Numerical methods
num_methods = {
    'spatialDerivs': 'SL6',  # SL4 or EX2
    'spatialDerivsBuffer': 'EX4',
    'timeStepping': 'RK4',  # RK4 or Euler
    'neumannOrder': 6,
    'neumann2Order': 2,
    'spatialFilterStrength': 0.49,  # -0.5 < alpha < 0.5
    'spatialFilterTime': 0,
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
# num_methods['SFD']['extraRegion'] = [{
#     'location': [(x1+x2)/2, 0, 0],
#     'size': [L, D, float('inf')],
#     'X': 0.1
# }]


