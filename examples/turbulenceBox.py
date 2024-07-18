import numpy as np

# Flow parameters
delta = 1
Re = 1000
Ma = 0.1
Lx = 1
Ly = Lx
Lz = Lx
nx = 64
ny = nx
nz = nx
epsilon = 1

# Case name
caseName = f"TurbulenceBox-Re{Re}-Ma{str(Ma).replace('.', '')}-eps{epsilon}-L{Lx}x{Ly}x{Lz}-n{nx}x{ny}x{nz}"

# Domain decomposition
p_row = 4
p_col = 4

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
    'xf': Lx,
    'yi': 0,
    'yf': Ly,
    'zi': 0,
    'zf': Lz
}

# Flow type
flowType = {
    'name': 'periodicBox',
    'initial': {
        'type': 'uniform',  # uniform, blasius or file
        'U0': 0,
        'addNoise': 1e-2
    }
}

e0 = 1/((flowParameters['gamma']**2 - flowParameters['gamma']) * flowParameters['Ma']**2)

flowType['disturb'] = [{
    'x': [-np.inf, np.inf],
    'y': [-np.inf, np.inf],
    'z': [-np.inf, np.inf],
    'var': 'UVWE',
    'type': 'turbulenceGenerator',
    'extraNodes': [0, 0, 0, 0, 0, 0],
    'par': [epsilon, e0],
    'active': True,
    'fitPoints': False,
    'forcing': True
}]

# Mesh parameters
mesh = {
    'x': {'n': nx, 'type': 'uniform', 'matchFixed': 2, 'periodic': True, 'fixPeriodicDomainSize': True, 'extraRefinement': 0},
    'y': {'n': ny, 'type': 'uniform', 'matchFixed': 2, 'periodic': True, 'fixPeriodicDomainSize': True, 'extraRefinement': 0},
    'z': {'n': nz, 'type': 'uniform', 'matchFixed': 2, 'periodic': True, 'fixPeriodicDomainSize': True, 'extraRefinement': 0}
}

for axis in ['x', 'y', 'z']:
    mesh[axis]['buffer'] = {'i': {'n': 0}, 'f': {'n': 0}}

# Tracked points
mesh['trackedPoints'] = [0, 0, 0]
mesh['trackedNorm'] = True

# Time control
time = {
    'control': 'cfl',
    'dt': 1,
    'maxCFL': 1.3,
    'qtimes': 1,
    'tmax': 100
}

logAll = 1

# Numerical methods
numMethods = {
    'spatialDerivs': 'SL6',  # SL4 or EX2
    'spatialDerivsBuffer': 'EX4',
    'timeStepping': 'RK4',  # RK4 or Euler
    'neumannOrder': 6,
    'neumann2Order': 2,
    'spatialFilterStrength': 0.499,  # -0.5 < alpha < 0.5
    'spatialFilterTime': 0,  # Characteristic time of the spatial filter (optional, default = 0)
    'filterDirections': [1, 1, 1],
    'filterBorders': 'reducedOrder',
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


