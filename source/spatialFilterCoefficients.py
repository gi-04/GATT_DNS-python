import numpy as np

def spatialFilterCoefficients(alpha, filterBorders):
    """
    This function outputs the coefficients for the 10th order spatial filter from Gaitonde 1998 for a given alpha.
    """
    
    filterStencilLHS = [1, alpha]
    
    filterStencilRHS = [
        (193 + 126 * alpha) / 256,
        (105 + 302 * alpha) / 512,
        (-15 + 30 * alpha) / 128,
        (45 - 90 * alpha) / 1024,
        (-5 + 10 * alpha) / 512,
        (1 - 2 * alpha) / 1024
    ]
    
    if isinstance(filterBorders, bool):
        filterBorders = 'decentered' if filterBorders else 'off'
    
    if filterBorders == 'decentered':
        filterDecenteredStencilLHS = [1, alpha]
        
        a_bound = np.zeros((11, 5))
        
        a_bound[0, 4] = (-1 + 2 * alpha) / 1024
        a_bound[1, 4] = (5 - 10 * alpha) / 512
        a_bound[2, 4] = (-45 + 90 * alpha) / 1024
        a_bound[3, 4] = (15 + 98 * alpha) / 128
        a_bound[4, 4] = (407 + 210 * alpha) / 512
        a_bound[5, 4] = (63 + 130 * alpha) / 256
        a_bound[6, 4] = (-105 + 210 * alpha) / 512
        a_bound[7, 4] = (15 - 30 * alpha) / 128
        a_bound[8, 4] = (-45 + 90 * alpha) / 1024
        a_bound[9, 4] = (5 - 10 * alpha) / 512
        a_bound[10, 4] = (-1 + 2 * alpha) / 1024

        a_bound[0, 3] = (1 - 2 * alpha) / 1024
        a_bound[1, 3] = (-5 + 10 * alpha) / 512
        a_bound[2, 3] = (45 + 934 * alpha) / 1024
        a_bound[3, 3] = (113 + 30 * alpha) / 128
        a_bound[4, 3] = (105 + 302 * alpha) / 512
        a_bound[5, 3] = (-63 + 126 * alpha) / 256
        a_bound[6, 3] = (105 - 210 * alpha) / 512
        a_bound[7, 3] = (-15 + 30 * alpha) / 128
        a_bound[8, 3] = (45 - 90 * alpha) / 1024
        a_bound[9, 3] = (-5 + 10 * alpha) / 512
        a_bound[10, 3] = (1 - 2 * alpha) / 1024

        a_bound[0, 2] = (-1 + 2 * alpha) / 1024
        a_bound[1, 2] = (5 + 502 * alpha) / 512
        a_bound[2, 2] = (979 + 90 * alpha) / 1024
        a_bound[3, 2] = (15 + 98 * alpha) / 128
        a_bound[4, 2] = (-105 + 210 * alpha) / 512
        a_bound[5, 2] = (63 - 126 * alpha) / 256
        a_bound[6, 2] = (-105 + 210 * alpha) / 512
        a_bound[7, 2] = (15 - 30 * alpha) / 128
        a_bound[8, 2] = (-45 + 90 * alpha) / 1024
        a_bound[9, 2] = (5 - 10 * alpha) / 512
        a_bound[10, 2] = (-1 + 2 * alpha) / 1024

        a_bound[0, 1] = (1 + 1022 * alpha) / 1024
        a_bound[1, 1] = (507 + 10 * alpha) / 512
        a_bound[2, 1] = (45 + 934 * alpha) / 1024
        a_bound[3, 1] = (-15 + 30 * alpha) / 128
        a_bound[4, 1] = (105 - 210 * alpha) / 512
        a_bound[5, 1] = (-63 + 126 * alpha) / 256
        a_bound[6, 1] = (105 - 210 * alpha) / 512
        a_bound[7, 1] = (-15 + 30 * alpha) / 128
        a_bound[8, 1] = (45 - 90 * alpha) / 1024
        a_bound[9, 1] = (-5 + 10 * alpha) / 512
        a_bound[10, 1] = (1 - 2 * alpha) / 1024

        a_bound[0, 0] = (1023 + 1 * alpha) / 1024
        a_bound[1, 0] = (5 + 507 * alpha) / 512
        a_bound[2, 0] = (-45 + 45 * alpha) / 1024
        a_bound[3, 0] = (15 - 15 * alpha) / 128
        a_bound[4, 0] = (-105 + 105 * alpha) / 512
        a_bound[5, 0] = (63 - 63 * alpha) / 256
        a_bound[6, 0] = (-105 + 105 * alpha) / 512
        a_bound[7, 0] = (15 - 15 * alpha) / 128
        a_bound[8, 0] = (-45 + 45 * alpha) / 1024
        a_bound[9, 0] = (5 - 5 * alpha) / 512
        a_bound[10, 0] = (-1 + 1 * alpha) / 1024
        
        filterDecenteredStencilRHS = a_bound.T
        
    elif filterBorders == 'reducedOrder':
        filterDecenteredStencilLHS = 1
        
        F2 = np.array([1/2 + alpha, 1/2 + alpha]) / np.array([1, 2])
        F4 = np.array([5/8 + 3/4 * alpha, 1/2 + alpha, -1/8 + 1/4 * alpha]) / np.array([1, 2, 2])
        F6 = np.array([11/16 + 5/8 * alpha, 15/32 + 17/16 * alpha, -3/16 + 3/8 * alpha, 1/32 - 1/16 * alpha]) / np.array([1, 2, 2, 2])
        F8 = np.array([93/128 + 70/128 * alpha, 7/16 + 18/16 * alpha, -7/32 + 14/32 * alpha, 1/16 - 1/8 * alpha, -1/128 + 1/64 * alpha]) / np.array([1, 2, 2, 2, 2])
        
        filterDecenteredStencilRHS = np.zeros((5, 9))
        
        filterDecenteredStencilRHS[0, 0] = 1
        filterDecenteredStencilRHS[1, 0:3] = F2[[1, 0, 1]]
        filterDecenteredStencilRHS[2, 0:5] = F4[[2, 1, 0, 1, 2]]
        filterDecenteredStencilRHS[3, 0:7] = F6[[3, 2, 1, 0, 1, 2, 3]]
        filterDecenteredStencilRHS[4, 0:9] = F8[[4, 3, 2, 1, 0, 1, 2, 3, 4]]
        
    elif filterBorders == 'off':
        filterDecenteredStencilLHS = np.diag(np.ones(5))
        filterDecenteredStencilRHS = np.diag(np.ones(5))
    
    return filterStencilLHS, filterStencilRHS, filterDecenteredStencilLHS, filterDecenteredStencilRHS



