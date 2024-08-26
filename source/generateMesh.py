import numpy as np
import scipy.interpolate as interp
import warnings

def generateMesh(xi, xf, mesh, direction):
    # This function computes the mesh for each direction and adds the buffer zones.

    # If the number of nodes is undefined, generate arbitrary mesh with same refining and define it
    if 'n' not in mesh:
        meshTemp = mesh.copy()
        meshTemp['n'] = int(np.ceil((xf - xi) / mesh['d0']))
        meshTemp['matchFixed'] = 0
        meshTemp['extraRefinement'] = 0
        meshTemp['fixPeriodicDomainSize'] = 0
        meshTemp['buffer']['i']['n'] = 0
        meshTemp['buffer']['f']['n'] = 0
        Xtemp, _, _ = generateMesh(xi, xf, meshTemp, direction)
        
        dBase = np.max(np.diff(Xtemp))
        
        mesh['n'] = int(np.ceil(meshTemp['n'] * dBase / mesh['d0']))

    # If the number of nodes in any buffer zone is undefined, define it based on the prescribed length
    if 'n' not in mesh['buffer']['i'] or 'n' not in mesh['buffer']['f']:
        meshTemp = mesh.copy()
        meshTemp['matchFixed'] = 0
        meshTemp['extraRefinement'] = 0
        meshTemp['fixPeriodicDomainSize'] = 0
        meshTemp['buffer']['i']['n'] = 0
        meshTemp['buffer']['f']['n'] = 0
        
        Xtemp, _, _ = generateMesh(xi, xf, meshTemp, direction)
        
        if 'n' not in mesh['buffer']['i']:
            dBase = Xtemp[1] - Xtemp[0]
            target = mesh['buffer']['i']['l']
            mesh['buffer']['i']['n'] = int(np.floor(target / dBase))
            XB = calcBufferZone(mesh['buffer']['i'])
            value = XB[-1] * dBase
            
            while value > target and mesh['buffer']['i']['n'] > 2:
                mesh['buffer']['i']['n'] = int(np.floor(mesh['buffer']['i']['n'] / 2))
                XB = calcBufferZone(mesh['buffer']['i'])
                value = XB[-1] * dBase
            while value < target:
                mesh['buffer']['i']['n'] += 1
                XB = calcBufferZone(mesh['buffer']['i'])
                value = XB[-1] * dBase
        
        if 'n' not in mesh['buffer']['f']:
            dBase = Xtemp[-1] - Xtemp[-2]
            target = mesh['buffer']['f']['l']
            mesh['buffer']['f']['n'] = int(np.floor(target / dBase))
            XB = calcBufferZone(mesh['buffer']['f'])
            value = XB[-1] * dBase
            
            while value > target and mesh['buffer']['f']['n'] > 2:
                mesh['buffer']['f']['n'] = int(np.floor(mesh['buffer']['f']['n'] / 2))
                XB = calcBufferZone(mesh['buffer']['f'])
                value = XB[-1] * dBase
            while value < target:
                mesh['buffer']['f']['n'] += 1
                XB = calcBufferZone(mesh['buffer']['f'])
                value = XB[-1] * dBase

    # If mesh is a single point, run as uniform
    if mesh['n'] == 1:
        mesh['type'] = 'uniform'

    # Remove buffer zone for periodic cases
    if mesh['periodic'] and (mesh['buffer']['i']['n'] > 0 or mesh['buffer']['f']['n'] > 0):
        warnings.warn(f'Buffer zone was removed for periodic dimension {direction}')
        mesh['buffer']['i']['n'] = 0
        mesh['buffer']['f']['n'] = 0

    # Add temporary node at the end if needed
    if mesh['periodic'] and mesh['fixPeriodicDomainSize'] and mesh['n'] > 1 and mesh['type'] != 'file':
        mesh['n'] += 1
        addedTempNode = True
    else:
        addedTempNode = False

    # Count nodes
    nx = mesh['n'] + mesh['buffer']['i']['n'] + mesh['buffer']['f']['n']
    physicalStart = mesh['buffer']['i']['n']
    physicalEnd = mesh['buffer']['i']['n'] + mesh['n']

    # Compute physical domain
    if mesh['type'] == 'uniform':
        XPhysical = np.linspace(xi, xf, mesh['n'])
    elif mesh['type'] == 'power':
        XPhysical = xi + (xf - xi) * np.linspace(0, 1, mesh['n']) ** mesh['power']
    elif mesh['type'] == 'tanh':
        if mesh['local'] == 'f':  # Refine at the end
            eta = np.tanh(np.linspace(0, mesh['par'], mesh['n']))
            eta = eta / eta[-1]
        elif mesh['local'] == 'i':  # Refine at the start
            eta = np.tanh(np.linspace(0, mesh['par'], mesh['n']))
            eta = eta / eta[-1]
            eta = 1 - eta[::-1]
        elif mesh['local'] == 'b':  # Refine at both sides
            eta = np.tanh(np.linspace(-mesh['par'], mesh['par'], mesh['n']))
            eta = eta / eta[-1]
            eta = (eta + 1) / 2
        XPhysical = (xf - xi) * eta + xi
    elif mesh['type'] == 'attractors':
        xBase = np.linspace(xi, xf, 100 * mesh['n'])
        eta = np.ones(100 * mesh['n'])
        for i in range(len(mesh['attractorPoints'])):
            eta += mesh['attractorStrength'][i] * np.exp(-((xBase - mesh['attractorPoints'][i]) / mesh['attractorSize'][i]) ** 2)

        # por ora decidi retirar todo este trecho porque, investigando as funções até aqui, não me parece que
        # o campo atractorRegions pode vir a existir em mesh, e sim apenas em mesh.x ou mesh.y - gigiaero (28/06/2024)
        # if 'attractorRegions' in mesh:
        #     for i in range(len(mesh['attractorRegions'])):            # troca do shape() por len() (gigiaero - 26/08/2024)
        #         nodePositions = mesh['attractorRegions'][i, :4]
        #         if np.isinf(nodePositions[1]):
        #             nodePositions[:2] = [xi - 2, xi - 1]
        #         if np.isinf(nodePositions[3]):
        #             nodePositions[2:] = [xf + 1, xf + 2]
        #         nodePositions = np.concatenate(([nodePositions[0] - 1], nodePositions, [nodePositions[3] + 1]))
        #         eta += mesh['attractorRegions'][i, 4] * interp.pchip_interpolate(nodePositions, [0, 0, 1, 1, 0, 0], xBase) ** 2


        eta = np.cumsum(eta)
        eta = eta - eta[0]
        eta = (mesh['n'] - 1) * eta / eta[-1] + 1
        eta[-1] = mesh['n']
        XPhysical = interp.pchip_interpolate(eta, xBase, np.arange(1, mesh['n'] + 1))
    elif mesh['type'] == 'attractors_old':
        xBase = np.linspace(xi, xf, mesh['n'])
        eta = np.ones(mesh['n'])
        for i in range(len(mesh['attractorPoints'])):
            eta += mesh['attractorStrength'][i] * np.exp(-((xBase - mesh['attractorPoints'][i]) / mesh['attractorSize'][i]) ** 2)
        eta = np.cumsum(eta)
        eta = eta - eta[0]
        eta = (mesh['n'] - 1) * eta / eta[-1] + 1
        eta[-1] = mesh['n']
        XPhysical = interp.pchip_interpolate(eta, xBase, np.arange(1, mesh['n'] + 1))
    elif mesh['type'] == 'file':
        X = np.loadtxt(mesh['file'])
        if mesh['file'].endswith('.mat'):
            X = X[direction]
        if X.shape[1] == 1:
            X = X.T
        if 'fileCalcBuffer' not in mesh or not mesh['fileCalcBuffer']:
            if len(X) != nx:
                raise ValueError(f'{direction} mesh from file {mesh["file"]} contains {len(X)} nodes instead of {nx} as specified in parameters')
            return X, mesh, nx
        if len(X) != mesh['n']:
            raise ValueError(f'{direction} mesh from file {mesh["file"]} contains {len(X)} nodes instead of {mesh["n"]} as specified in parameters')
        XPhysical = X

    # Match fixed points such as cavity edges
    if mesh['matchFixed'] and mesh['n'] > 1:
        fixPoints = np.unique(np.concatenate(([XPhysical[0], XPhysical[1]], mesh['fixPoints'], [XPhysical[-2], XPhysical[-1]])))
        closestNodes = np.argmin(np.abs(XPhysical[:, None] - fixPoints), axis=0)
        ok = not np.any(np.diff(closestNodes) == 0)  # Check for two fixed points assigned to the same node
        iter = 0
        iterMax = mesh['n']
        while not ok and iter <= iterMax:
            duplicates = np.where(np.diff(closestNodes) == 0)[0]
            closestNodes[duplicates] -= 1
            ok = not np.any(np.diff(closestNodes) == 0)
            iter += 1
        if mesh['matchFixed'] == 2:
            XPhysical -= interp.pchip_interpolate(closestNodes, XPhysical[closestNodes] - fixPoints, np.arange(mesh['n']))
        else:
            XPhysical -= interp.spline_interpolate(closestNodes, XPhysical[closestNodes] - fixPoints, np.arange(mesh['n']))
        XPhysical[closestNodes] = fixPoints  # Done to fix any possible rounding errors

    X = np.zeros(nx)
    X[physicalStart:physicalEnd] = XPhysical

    # Add buffer zones
    if mesh['buffer']['i']['n'] > 0:
        baseDist = X[physicalStart + 1] - X[physicalStart]
        XB = calcBufferZone(mesh['buffer']['i']) * baseDist
        X[:physicalStart] = X[physicalStart] - XB[::-1]

    if mesh['buffer']['f']['n'] > 0:
        baseDist = X[physicalEnd - 1] - X[physicalEnd - 2]
        XB = calcBufferZone(mesh['buffer']['f']) * baseDist
        X[physicalEnd:] = X[physicalEnd - 1] + XB

    # Add extra refinement nodes
    if mesh['extraRefinement'] > 0:
        # Recalculate number of nodes
        er = mesh['extraRefinement']
        mesh['n'] = (1 + er) * mesh['n'] - er
        mesh['buffer']['i']['n'] = (1 + er) * mesh['buffer']['i']['n']
        mesh['buffer']['f']['n'] = (1 + er) * mesh['buffer']['f']['n']
        nx = mesh['n'] + mesh['buffer']['i']['n'] + mesh['buffer']['f']['n']
        
        # Add new nodes
        Xnew = interp.pchip_interpolate(np.linspace(0, 1, len(X)), X, np.linspace(0, 1, nx))
        Xnew[::(1 + er)] = X
        X = Xnew

    # Remove temporary node
    if addedTempNode:
        X = X[:-1]
        mesh['n'] -= 1
        nx -= 1

    return X, mesh, nx

def calcBufferZone(par):
    if par['type'] == 'exponential':
        if 'ramp' not in par:
            XB = [1 + par['stretching']]
            for i in range(1, par['n']):
                XB.append(XB[-1] + (1 + par['stretching']) ** i)
        else:
            delta = np.ones(par['n'])
            for i in range(1, par['n']):
                stretching = 1 + min(1, i / par['ramp']) * par['stretching']
                delta[i] = delta[i - 1] * stretching
            XB = np.cumsum(delta)
    elif par['type'] == 'sigmoid':
        lambda_ = 1
        xb = np.linspace(-10, 2, par['n'])
        delta = 1 / (1 + np.exp(-lambda_ * xb))
        delta = delta * par['stretching'] + 1
        XB = np.cumsum(delta)
    return XB



