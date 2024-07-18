
import numpy as np

def init_boundaries(boundary, mesh, domain_slices_y, domain_slices_z, p_row, p_col):
    # This function initializes the boundary conditions to the domain
    # Before running the DNS, it should be called with the boundary structure
    # as input

    # Variable names are as follows:
    # First letter - What it represents
    #   n is the number of boundaries of that type
    #   i is the index this boundary is located
    #   v is the value a Dirichlet boundary will be set to
    #   d is the direction a null Neumann boundary will be applied
    # Second letter - The flow variable it applies to
    # Third letter - The type of boundary: d for Dirichlet and n for Neumann

    # One output is created for each domain slice

    # Loop through all boundaries and organize them in the correct variables
    biG = {
        'nUd': 0, 'nVd': 0, 'nWd': 0, 'nPd': 0, 'nEd': 0,
        'nUn': 0, 'nVn': 0, 'nWn': 0, 'nPn': 0, 'nEn': 0,
        'nUs': 0, 'nVs': 0, 'nWs': 0, 'nPs': 0, 'nEs': 0,
        'iUd': np.empty((0, 6), dtype=int), 'iVd': np.empty((0, 6), dtype=int), 'iWd': np.empty((0, 6), dtype=int), 'iPd': np.empty((0, 6), dtype=int), 'iEd': np.empty((0, 6), dtype=int),
        'iUn': np.empty((0, 6), dtype=int), 'iVn': np.empty((0, 6), dtype=int), 'iWn': np.empty((0, 6), dtype=int), 'iPn': np.empty((0, 6), dtype=int), 'iEn': np.empty((0, 6), dtype=int),
        'iUs': np.empty((0, 6), dtype=int), 'iVs': np.empty((0, 6), dtype=int), 'iWs': np.empty((0, 6), dtype=int), 'iPs': np.empty((0, 6), dtype=int), 'iEs': np.empty((0, 6), dtype=int),
        'vUd': np.empty(0), 'vVd': np.empty(0), 'vWd': np.empty(0), 'vPd': np.empty(0), 'vEd': np.empty(0),
        'dUn': np.empty(0), 'dVn': np.empty(0), 'dWn': np.empty(0), 'dPn': np.empty(0), 'dEn': np.empty(0),
        'dUs': np.empty(0), 'dVs': np.empty(0), 'dWs': np.empty(0), 'dPs': np.empty(0), 'dEs': np.empty(0),
    }

    direction_order = ['xi', 'xf', 'yi', 'yf', 'zi', 'zf']
    for i in range(len(boundary['val'])):
        if boundary['type'][i] == 'dir':
            if boundary['var'][i] == 'u':
                biG['nUd'] += 1
                biG['iUd'] = np.vstack((biG['iUd'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['vUd'] = np.append(biG['vUd'], boundary['val'][i])
            elif boundary['var'][i] == 'v':
                biG['nVd'] += 1
                biG['iVd'] = np.vstack((biG['iVd'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['vVd'] = np.append(biG['vVd'], boundary['val'][i])
            elif boundary['var'][i] == 'w':
                biG['nWd'] += 1
                biG['iWd'] = np.vstack((biG['iWd'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['vWd'] = np.append(biG['vWd'], boundary['val'][i])
            elif boundary['var'][i] == 'p':
                biG['nPd'] += 1
                biG['iPd'] = np.vstack((biG['iPd'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['vPd'] = np.append(biG['vPd'], boundary['val'][i])
            elif boundary['var'][i] == 'e':
                biG['nEd'] += 1
                biG['iEd'] = np.vstack((biG['iEd'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['vEd'] = np.append(biG['vEd'], boundary['val'][i])
        elif boundary['type'][i] == 'neu':
            if boundary['var'][i] == 'u':
                biG['nUn'] += 1
                biG['iUn'] = np.vstack((biG['iUn'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['dUn'] = np.append(biG['dUn'], np.where(np.array(direction_order) == boundary['dir'][i])[0][0])
            elif boundary['var'][i] == 'v':
                biG['nVn'] += 1
                biG['iVn'] = np.vstack((biG['iVn'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['dVn'] = np.append(biG['dVn'], np.where(np.array(direction_order) == boundary['dir'][i])[0][0])
            elif boundary['var'][i] == 'w':
                biG['nWn'] += 1
                biG['iWn'] = np.vstack((biG['iWn'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['dWn'] = np.append(biG['dWn'], np.where(np.array(direction_order) == boundary['dir'][i])[0][0])
            elif boundary['var'][i] == 'p':
                biG['nPn'] += 1
                biG['iPn'] = np.vstack((biG['iPn'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['dPn'] = np.append(biG['dPn'], np.where(np.array(direction_order) == boundary['dir'][i])[0][0])
            elif boundary['var'][i] == 'e':
                biG['nEn'] += 1
                biG['iEn'] = np.vstack((biG['iEn'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['dEn'] = np.append(biG['dEn'], np.where(np.array(direction_order) == boundary['dir'][i])[0][0])
        elif boundary['type'][i] == 'sec':
            if boundary['var'][i] == 'u':
                biG['nUs'] += 1
                biG['iUs'] = np.vstack((biG['iUs'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['dUs'] = np.append(biG['dUs'], np.where(np.array(direction_order) == boundary['dir'][i])[0][0])
            elif boundary['var'][i] == 'v':
                biG['nVs'] += 1
                biG['iVs'] = np.vstack((biG['iVs'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['dVs'] = np.append(biG['dVs'], np.where(np.array(direction_order) == boundary['dir'][i])[0][0])
            elif boundary['var'][i] == 'w':
                biG['nWs'] += 1
                biG['iWs'] = np.vstack((biG['iWs'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['dWs'] = np.append(biG['dWs'], np.where(np.array(direction_order) == boundary['dir'][i])[0][0])
            elif boundary['var'][i] == 'p':
                biG['nPs'] += 1
                biG['iPs'] = np.vstack((biG['iPs'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['dPs'] = np.append(biG['dPs'], np.where(np.array(direction_order) == boundary['dir'][i])[0][0])
            elif boundary['var'][i] == 'e':
                biG['nEs'] += 1
                biG['iEs'] = np.vstack((biG['iEs'], np.array([boundary['xi'][i], boundary['xf'][i], boundary['yi'][i], boundary['yf'][i], boundary['zi'][i], boundary['zf'][i]])))
                biG['dEs'] = np.append(biG['dEs'], np.where(np.array(direction_order) == boundary['dir'][i])[0][0])

    # Get information on corners
    biG['cL'] = boundary['corners']['limits']
    biG['cD'] = boundary['corners']['dir']
    biG['adiabatic'] = boundary['corners']['adiabatic']
    biG['cN'] = len(biG['cL'])

    neumann_length = len(boundary['neumannCoeffs'])
    neumann2_length = len(boundary['neumann2Coeffs'])

    # Set gamma-1 value for converting pressure to density
    biG['gamma1'] = boundary['gamma'] - 1

    # Find regions inside walls
    #biG['wR'] = ~boundary['flowRegion']
    biG['E0'] = boundary['E0']

    # Split boundaries for different processors
    bi = [None] * (p_row * p_col)
    for j in range(p_row):
        for k in range(p_col):
            n_proc = k + (j - 1) * p_col

            biL = biG.copy()

            Ji = domain_slices_y[0, j]
            Jf = domain_slices_y[1, j]
            Ki = domain_slices_z[0, k]
            Kf = domain_slices_z[1, k]

            biL['iUd'], biL['nUd'], biL['vUd'] = limit_indices(biL['iUd'], biL['nUd'], biL['vUd'], 'd', Ji, Jf, Ki, Kf, 0)
            biL['iVd'], biL['nVd'], biL['vVd'] = limit_indices(biL['iVd'], biL['nVd'], biL['vVd'], 'd', Ji, Jf, Ki, Kf, 0)
            biL['iWd'], biL['nWd'], biL['vWd'] = limit_indices(biL['iWd'], biL['nWd'], biL['vWd'], 'd', Ji, Jf, Ki, Kf, 0)
            biL['iPd'], biL['nPd'], biL['vPd'] = limit_indices(biL['iPd'], biL['nPd'], biL['vPd'], 'd', Ji, Jf, Ki, Kf, 0)
            biL['iEd'], biL['nEd'], biL['vEd'] = limit_indices(biL['iEd'], biL['nEd'], biL['vEd'], 'd', Ji, Jf, Ki, Kf, 0)

            biL['iUn'], biL['nUn'], biL['dUn'] = limit_indices(biL['iUn'], biL['nUn'], biL['dUn'], 'n', Ji, Jf, Ki, Kf, neumann_length)
            biL['iVn'], biL['nVn'], biL['dVn'] = limit_indices(biL['iVn'], biL['nVn'], biL['dVn'], 'n', Ji, Jf, Ki, Kf, neumann_length)
            biL['iWn'], biL['nWn'], biL['dWn'] = limit_indices(biL['iWn'], biL['nWn'], biL['dWn'], 'n', Ji, Jf, Ki, Kf, neumann_length)
            biL['iPn'], biL['nPn'], biL['dPn'] = limit_indices(biL['iPn'], biL['nPn'], biL['dPn'], 'n', Ji, Jf, Ki, Kf, neumann_length)
            biL['iEn'], biL['nEn'], biL['dEn'] = limit_indices(biL['iEn'], biL['nEn'], biL['dEn'], 'n', Ji, Jf, Ki, Kf, neumann_length)

            biL['iUs'], biL['nUs'], biL['dUs'] = limit_indices(biL['iUs'], biL['nUs'], biL['dUs'], 'n', Ji, Jf, Ki, Kf, neumann2_length)
            biL['iVs'], biL['nVs'], biL['dVs'] = limit_indices(biL['iVs'], biL['nVs'], biL['dVs'], 'n', Ji, Jf, Ki, Kf, neumann2_length)
            biL['iWs'], biL['nWs'], biL['dWs'] = limit_indices(biL['iWs'], biL['nWs'], biL['dWs'], 'n', Ji, Jf, Ki, Kf, neumann2_length)
            biL['iPs'], biL['nPs'], biL['dPs'] = limit_indices(biL['iPs'], biL['nPs'], biL['dPs'], 'n', Ji, Jf, Ki, Kf, neumann2_length)
            biL['iEs'], biL['nEs'], biL['dEs'] = limit_indices(biL['iEs'], biL['nEs'], biL['dEs'], 'n', Ji, Jf, Ki, Kf, neumann2_length)

            values = np.concatenate((biL['cD'], biL['adiabatic'].reshape((-1, 1))), axis=1)
            biL['cL'], biL['cN'], values = limit_indices(biL['cL'], biL['cN'], values, 'c', Ji, Jf, Ki, Kf, neumann_length)
            biL['cD'] = values[:, :3]
            biL['adiabatic'] = values[:, 3]

            bi[n_proc] = biL

    # Split disturbances across processors
    for j in range(p_row):
        for k in range(p_col):
            n_proc = k + (j - 1) * p_col
            disturb = []
            for i in range(len(boundary['disturb'])):
                if boundary['disturb'][i] is not None:
                    ind = boundary['disturb'][i]['ind']
                    ind[2:5] = [max(ind[2], domain_slices_y[0, j]), min(ind[3], domain_slices_y[1, j]), max(ind[4], domain_slices_z[0, k]), min(ind[5], domain_slices_z[1, k])]

                    if ind[2] <= ind[3] and ind[4] <= ind[5]:
                        disturb.append(boundary['disturb'][i])
                        disturb[-1]['ind'] = ind
                        disturb[-1]['X'] = mesh['X'][ind[0]:ind[1]]
                        disturb[-1]['Y'] = mesh['Y'][ind[2]:ind[3]]
                        disturb[-1]['Z'] = mesh['Z'][ind[4]:ind[5]]

            bi[n_proc]['disturb'] = disturb

    return bi

def limit_indices(ind, n, vd, type, Ji, Jf, Ki, Kf, neumann_length):
    if n == 0:
        return ind, n, vd

    ind[:, [2, 3, 4, 5]] = np.maximum(ind[:, [2, 3, 4, 5]], np.array([Ji, Ki, Ji, Ki]))
    ind[:, [2, 3, 4, 5]] = np.minimum(ind[:, [2, 3, 4, 5]], np.array([Jf, Kf, Jf, Kf]))

    to_remove = np.any(ind[:, [2, 3, 4, 5]] > ind[:, [3, 4, 5, 6]], axis=1)
    ind = np.delete(ind, to_remove, axis=0)

    if type == 'n':
        if len(vd) == 1 and n > 1:
            vd = np.delete(vd, to_remove)
        else:
            vd = np.delete(vd, to_remove, axis=0)

        for i in range(n - sum(to_remove)):
            switch = vd[i]
            if switch == 3:
                if ind[i, 3] + neumann_length > Jf:
                    raise ValueError(f'There is a y+ Neumann condition at J = {ind[i, 3]} crossing a domain slice at J = {Jf}. Consider changing p_row.')
            elif switch == 4:
                if ind[i, 2] - neumann_length < Ji:
                    raise ValueError(f'There is a y- Neumann condition at J = {ind[i, 2]} crossing a domain slice at J = {Ji}. Consider changing p_row.')
            elif switch == 5:
                if ind[i, 5] + neumann_length > Kf:
                    raise ValueError(f'There is a z+ Neumann condition at K = {ind[i, 5]} crossing a domain slice at K = {Kf}. Consider changing p_col.')
            elif switch == 6:
                if ind[i, 4] - neumann_length < Ki:
                    raise ValueError(f'There is a z- Neumann condition at K = {ind[i, 4]} crossing a domain slice at K = {Ki}. Consider changing p_col.')

    return ind, n - sum(to_remove), vd


