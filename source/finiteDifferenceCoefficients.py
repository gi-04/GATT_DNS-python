import numpy as np

def finiteDifferenceCoefficients(method):
    """
    This function contains the coefficients for each finite differences method. 
    New methods can be defined here if needed.

    Centered Stencils start at the center point and move to the right hand side.
    Uncentered stencils are as they should be in the left hand side top of the matrix.
    """
    if method == 'SL6':
        w = 1.8

        # Center
        matriz_a = np.array([
            [1, 1, 1, -2],
            [1, 4, 9, -6],
            [1, 16, 81, -10],
            [np.sin(w), np.sin(2*w)/2, np.sin(3*w)/3, -2*w*np.cos(w)]
        ])
        matriz_b = np.array([1, 0, 0, w])
        coeffs = np.linalg.solve(matriz_a, matriz_b)

        a = coeffs[0] / 2
        b = coeffs[1] / 4
        c = coeffs[2] / 6
        alpha = coeffs[3]

        centeredStencilLHS = [1, alpha]
        centeredStencilRHS = [0, a, b, c]

        # Border
        # First stage
        matriz_a = np.array([
            [1, 1, 1, 1, 1, 1, 0],
            [0, 1, 2, 3, 4, 5, -1],
            [0, 1, 4, 9, 16, 25, -2],
            [0, 1, 8, 27, 64, 125, -3],
            [0, 1, 16, 81, 256, 625, -4],
            [0, 1, 32, 243, 1024, 3125, -5],
            [0, 1, 64, 729, 4096, 15625, -6]
        ])
        matriz_b = np.array([0, 1, 0, 0, 0, 0, 0])
        coeffs = np.linalg.solve(matriz_a, matriz_b)

        a = coeffs[0]
        b = coeffs[1]
        c = coeffs[2]
        d = coeffs[3]
        e = coeffs[4]
        f = coeffs[5]
        alpha = coeffs[6]

        matriz_lhs_aux = [1, alpha, 0, 0]
        matriz_rhs_aux = [a, b, c, d, e, f]

        # Second stage
        matriz_a = np.array([
            [-1, 2],
            [-1, 6]
        ])
        matriz_b = np.array([-1, 0])
        coeffs = np.linalg.solve(matriz_a, matriz_b)

        a = coeffs[0] / 2
        alpha = coeffs[1]

        matriz_lhs_aux2 = [alpha, 1, alpha, 0]
        matriz_rhs_aux2 = [-a, 0, a, 0, 0, 0]

        # Third stage (same as centered SL4)
        matriz_a = np.array([
            [1, 0, -2/3],
            [0, 1, -4/3],
            [np.sin(w), np.sin(2*w)/2, -2*w*np.cos(w)]
        ])
        matriz_b = np.array([4/3, -1/3, w])
        coeffs = np.linalg.solve(matriz_a, matriz_b)

        a = coeffs[0] / 2
        b = coeffs[1] / 4
        alpha = coeffs[2]

        matriz_lhs_aux3 = [0, alpha, 1, alpha]
        matriz_rhs_aux3 = [-b, -a, 0, a, b, 0]

        decenteredStencilLHS = np.array([matriz_lhs_aux, matriz_lhs_aux2, matriz_lhs_aux3])
        decenteredStencilRHS = np.array([matriz_rhs_aux, matriz_rhs_aux2, matriz_rhs_aux3])

    elif method == 'SL6O3':
        # Inner-Scheme: SL6 - Optimized, alpha determined minimizing the integral of dispersion
        # from w = [0,wf], i.e., alpha = f(integral(Edisp)|_0^wf)), wf = pi/2;
        # P3: SL4 - Optimized following the procedure above, wf = pi/2;
        # P2: C4 - Pade Scheme
        # P1: C3 - Optimized, max. dissipation error Ediss(w=pi) = pi/2; (Edisp(w=pi/2)~6%)
        coeffs = [0.392465753424658, 1.565410958904110, 0.237260273972603, -0.017739726027397]
        alpha = coeffs[0]
        a = coeffs[1]
        b = coeffs[2]
        c = coeffs[3]
        centeredStencilLHS = [1, alpha]
        centeredStencilRHS = [0, a/2, b/4, c/6]

        coeffs_P3 = [0.350978473581213, 1.567318982387476, 0.134637964774951]
        alpha_P3 = coeffs_P3[0]
        a_P3 = coeffs_P3[1]
        b_P3 = coeffs_P3[2]

        decenteredStencilLHS = np.array([
            [1, (3*np.pi + 40)/(3*np.pi + 8), 0, 0],
            [1/4, 1, 1/4, 0],
            [0, alpha_P3, 1, alpha_P3]
        ])
        decenteredStencilRHS = np.array([
            [-(13*np.pi + 56)/(2*(3*np.pi + 8)), (15*np.pi + 8)/(2*(3*np.pi + 8)), -(3*np.pi - 56)/(2*(3*np.pi + 8)), (np.pi - 8)/(2*(3*np.pi + 8)), 0],
            [-3/4, 0, 3/4, 0, 0],
            [-b_P3/4, -a_P3/2, 0, a_P3/2, b_P3/4]
        ])

    elif method == 'SL4':
        # 4th order spectral-like compact derivatives
        w = 1.8

        matriz_a = np.array([
            [1, 0, -2/3],
            [0, 1, -4/3],
            [np.sin(w), np.sin(2*w)/2, -2*w*np.cos(w)]
        ])
        matriz_b = np.array([4/3, -1/3, w])
        coeffs = np.linalg.solve(matriz_a, matriz_b)

        a = coeffs[0] / 2
        b = coeffs[1] / 4
        alpha = coeffs[2]

        centeredStencilLHS = [1, alpha]
        centeredStencilRHS = [0, a, b]

        decenteredStencilLHS = np.array([
            [1, 3, 0],
            [1/4, 1, 1/4]
        ])
        decenteredStencilRHS = np.array([
            [-17/6, 3/2, 3/2, -1/6],
            [-3/4, 0, 3/4, 0]
        ])

    elif method == 'EX2':
        # 2nd order explicit finite differences
        centeredStencilLHS = 1
        centeredStencilRHS = [0, 1/2]

        decenteredStencilLHS = 1
        decenteredStencilRHS = [-3/2, 2, -1/2]

    elif method == 'EX4':
        # 4th order explicit finite differences
        centeredStencilLHS = 1
        centeredStencilRHS = [0, 2/3, -1/12]

        decenteredStencilLHS = 1
        decenteredStencilRHS = np.array([
            [-25/12, 4, -3, 4/3, -1/4],
            [-1/2, 0, 1/2, 0, 0]
        ])

    else:
        raise ValueError(f'Finite differences method not implemented: {method}. Check finiteDifferenceCoefficients file for available methods')

    return centeredStencilLHS, centeredStencilRHS, decenteredStencilLHS, decenteredStencilRHS



