import numpy as np

def calcCompressibleBL(flowParameters, adiabWall, mesh):
    Re = flowParameters['Re']  # Reynolds Number
    Minf = flowParameters['Ma']  # Mach Number
    Pr = flowParameters['Pr']  # Prandtl Number
    Tinf = flowParameters['T0']  # Temp. [K]
    gamma = flowParameters['gamma']  # Gamma

    xR = Re / 1.7208**2  # Reference x (where delta*(xR) = 1, assuming incompressible Blasius)
    E0 = 1 / (gamma * (gamma - 1) * Minf**2)

    Twall = 1  # Twall/Tinf
    if 'Twall' in flowParameters:
        Twall = flowParameters['Twall'] / Tinf  # Twall/Tinf
    sol = solve_compressibleBL(Minf, Pr, Tinf, Twall, gamma, adiabWall)

    eta_xR = sol['eta']
    U_xR = sol['f_p']
    R_xR = 1 / sol['rbar']
    E_xR = sol['rbar'] * E0
    V_xR = -(sol['rbar']) * (sol['f'] - eta_xR * sol['f_p']) * (1.7208 / np.sqrt(2)) * (1 / Re)

    y_xR = eta_xR * np.sqrt(2) / 1.7208
    dS_xR = np.trapz(y_xR, 1 - U_xR * R_xR / (U_xR[-1] * R_xR[-1]))

    print(' Initial Flow: Compressible BL')
    if adiabWall:
        print('   Adiabatic Wall, (dT/dy)_wall = 0, (Twall/TInf) = %1.3f' % (E_xR[0] / E0))
    else:
        print('   Isothermal Wall, T_wall = cte,    (Twall/TInf) = %1.3f' % (E_xR[0] / E0))

    print('   !!! BL thickness, (effective) / (incomp. BL, Blasius) = %1.4f' % dS_xR)

    # --------------------------
    X = mesh['X']
    nx = len(X)
    Y = mesh['Y']
    ny = len(Y)
    Z = mesh['Z']
    nz = len(Z)

    R = np.ones((nx, ny))
    U = np.ones((nx, ny))
    V = np.zeros((nx, ny))
    W = np.zeros((nx, ny))
    E = np.ones((nx, ny)) * E0

    indX = np.where(X > 0)[0]
    indY = np.where(Y >= 0)[0]
    for ix in range(len(indX)):
        y_xL = y_xR * np.sqrt(X[indX[ix]] / xR)

        U[indX[ix], indY] = np.interp(Y[indY], y_xL, U_xR, left=U_xR[-1], right=U_xR[-1])
        V[indX[ix], indY] = np.interp(Y[indY], y_xL, V_xR, left=V_xR[-1], right=V_xR[-1])
        R[indX[ix], indY] = np.interp(Y[indY], y_xL, R_xR, left=R_xR[-1], right=R_xR[-1])
        E[indX[ix], indY] = np.interp(Y[indY], y_xL, E_xR, left=E_xR[-1], right=E_xR[-1])

    U[X == 0, Y == 0] = 0
    indY = np.where(Y < 0)[0]
    U[:, indY] = np.tile(U[:, Y == 0], (1, len(indY)))
    V[:, indY] = np.tile(V[:, Y == 0], (1, len(indY)))
    R[:, indY] = np.tile(R[:, Y == 0], (1, len(indY)))
    E[:, indY] = np.tile(E[:, Y == 0], (1, len(indY)))

    U = np.repeat(U[:, :, np.newaxis], nz, axis=2)
    V = np.repeat(V[:, :, np.newaxis], nz, axis=2)
    W = np.repeat(W[:, :, np.newaxis], nz, axis=2)
    R = np.repeat(R[:, :, np.newaxis], nz, axis=2)
    E = np.repeat(E[:, :, np.newaxis], nz, axis=2)

    flow = {'U': U, 'V': V, 'W': W, 'R': R, 'E': E}
    return flow

def solve_compressibleBL(Minf, Pr, Tinf, Twall, Gamma, adiabWall):
    C2 = 110  # Sutherland Coefficient [Kelvin]
    lim = 10  # The value which simulates lim -> inf
    N = 500  # Number of Point
    h = lim / N  # Delta y
    delta = 1e-10  # Small Number for shooting method
    eps = 1e-9

    adi = 1 if adiabWall else 0

    # Initializing
    y1 = np.zeros(N + 1)  # f
    y2 = np.zeros(N + 1)  # f'
    y3 = np.zeros(N + 1)  # f''
    y4 = np.zeros(N + 1)  # rho(eta)
    y5 = np.zeros(N + 1)  # rho(eta)'
    eta = np.linspace(0, lim, N + 1)  # Iteration of eta up to infinity
    dalfa = 0
    dbeta = 0

    if adi == 1:
        # Boundary Conditions for Adiabatic Case
        y1[0] = 0
        y2[0] = 0
        y5[0] = 0

        # Initial Guess for the beginning of simulation
        alfa0 = 0.1  # Initial Guess
        beta0 = 3  # Initial Guess
    else:
        # Boundary Conditions for Isothermal Case
        y1[0] = 0
        y2[0] = 0
        y4[0] = Twall

        # Initial Guess for Beginning of Simulation
        alfa0 = 0.1  # Initial Guess
        beta0 = 3  # Initial Guess

    for ite in range(100000):
        if adi == 1:
            # Boundary Conditions for Adiabatic Case
            y1[0] = 0
            y2[0] = 0
            y5[0] = 0

            y3[0] = alfa0
            y4[0] = beta0
        else:
            # Boundary Conditions for Isothermal Case
            y1[0] = 0
            y2[0] = 0
            y4[0] = Twall

            y3[0] = alfa0
            y5[0] = beta0

        y1, y2, y3, y4, y5 = RK(eta, h, y1, y2, y3, y4, y5, C2, Tinf, Minf, Pr, Gamma)

        y2old = y2[-1]
        y4old = y4[-1]

        if adi == 1:
            # Boundary Conditions for Adiabatic Case
            y1[0] = 0
            y2[0] = 0
            y5[0] = 0

            y3[0] = alfa0 + delta
            y4[0] = beta0
        else:
            # Boundary Conditions for Isothermal Case
            y1[0] = 0
            y2[0] = 0
            y4[0] = Twall

            y3[0] = alfa0 + delta
            y5[0] = beta0

        y1, y2, y3, y4, y5 = RK(eta, h, y1, y2, y3, y4, y5, C2, Tinf, Minf, Pr, Gamma)

        y2new1 = y2[-1]
        y4new1 = y4[-1]

        if adi == 1:
            # Boundary Conditions for Adiabatic Case
            y1[0] = 0
            y2[0] = 0
            y5[0] = 0

            y3[0] = alfa0
            y4[0] = beta0 + delta
        else:
            # Boundary Conditions for Isothermal Case
            y1[0] = 0
            y2[0] = 0
            y4[0] = Twall

            y3[0] = alfa0
            y5[0] = beta0 + delta

        y1, y2, y3, y4, y5 = RK(eta, h, y1, y2, y3, y4, y5, C2, Tinf, Minf, Pr, Gamma)

        y2new2 = y2[-1]
        y4new2 = y4[-1]

        a11 = (y2new1 - y2old) / delta
        a21 = (y4new1 - y4old) / delta
        a12 = (y2new2 - y2old) / delta
        a22 = (y4new2 - y4old) / delta
        r1 = 1 - y2old
        r2 = 1 - y4old
        dalfa = (a22 * r1 - a12 * r2) / (a11 * a22 - a12 * a21)
        dbeta = (a11 * r2 - a21 * r1) / (a11 * a22 - a12 * a21)
        alfa0 += dalfa
        beta0 += dbeta

        if (abs(y2[-1] - 1) < eps) and (abs(y4[-1] - 1) < eps):
            Truey2 = y2[0]
            Truey4 = y4[0]
            break

    xaxis = np.zeros(len(eta))
    for i in range(1, len(eta)):
        xaxis[i] = (eta[i] - 0) * (y4[0] + 2 * np.sum(y4[1:i]) + y4[i]) / (2 * eta[i]) * h

    sol = {
        'eta': xaxis,
        'f': y1,
        'f_p': y2,
        'f_pp': y3,
        'rbar': y4,
        'rbar_p': y5
    }
    return sol

def Y1(y2):
    return

def Y2(y3):
    return

def Y3(y1, y3, y4, y5, C2, Tinf):
    RHS = -y3 * ((y5 / (2 * (y4))) - (y5 / (y4 + C2 / Tinf))) \
          - y1 * y3 * ((y4 + C2 / Tinf) / (np.sqrt(y4) * (1 + C2 / Tinf)))
    return RHS

def Y4(y5):
    return

def Y5(y1, y3, y4, y5, C2, Tinf, Minf, Pr, Gamma):
    RHS = -y5**2 * ((0.5 / y4) - (1 / (y4 + C2 / Tinf))) \
          - Pr * y1 * y5 / np.sqrt(y4) * (y4 + C2 / Tinf) / (1 + C2 / Tinf) \
          - (Gamma - 1) * Pr * Minf**2 * y3**2
    return RHS

def RK(eta, h, y1, y2, y3, y4, y5, C2, Tinf, Minf, Pr, Gamma):
    for i in range(len(eta) - 1):
        k11 = Y1(y2[i])
        k21 = Y2(y3[i])
        k31 = Y3(y1[i], y3[i], y4[i], y5[i], C2, Tinf)
        k41 = Y4(y5[i])
        k51 = Y5(y1[i], y3[i], y4[i], y5[i], C2, Tinf, Minf, Pr, Gamma)

        k12 = Y1(y2[i] + 0.5 * h * k21)
        k22 = Y2(y3[i] + 0.5 * h * k31)
        k32 = Y3(y1[i] + 0.5 * h * k11, y3[i] + 0.5 * h * k31, y4[i] + 0.5 * h * k41, y5[i] + 0.5 * h * k51, C2, Tinf)
        k42 = Y4(y5[i] + 0.5 * h * k51)
        k52 = Y5(y1[i] + 0.5 * h * k11, y3[i] + 0.5 * h * k31, y4[i] + 0.5 * h * k41, y5[i] + 0.5 * h * k51, C2, Tinf, Minf, Pr, Gamma)

        k13 = Y1(y2[i] + 0.5 * h * k22)
        k23 = Y2(y3[i] + 0.5 * h * k32)
        k33 = Y3(y1[i] + 0.5 * h * k12, y3[i] + 0.5 * h * k32, y4[i] + 0.5 * h * k42, y5[i] + 0.5 * h * k52, C2, Tinf)
        k43 = Y4(y5[i] + 0.5 * h * k52)
        k53 = Y5(y1[i] + 0.5 * h * k12, y3[i] + 0.5 * h * k32, y4[i] + 0.5 * h * k42, y5[i] + 0.5 * h * k52, C2, Tinf, Minf, Pr, Gamma)

        k14 = Y1(y2[i] + h * k23)
        k24 = Y2(y3[i] + h * k33)
        k34 = Y3(y1[i] + h * k13, y3[i] + h * k33, y4[i] + h * k43, y5[i] + h * k53, C2, Tinf)
        k44 = Y4(y5[i] + h * k53)
        k54 = Y5(y1[i] + h * k13, y3[i] + h * k33, y4[i] + h * k43, y5[i] + h * k53, C2, Tinf, Minf, Pr, Gamma)

        y5[i + 1] = y5[i] + (1 / 6) * (k51 + 2 * k52 + 2 * k53 + k54) * h
        y4[i + 1] = y4[i] + (1 / 6) * (k41 + 2 * k42 + 2 * k43 + k44) * h
        y3[i + 1] = y3[i] + (1 / 6) * (k31 + 2 * k32 + 2 * k33 + k34) * h
        y2[i + 1] = y2[i] + (1 / 6) * (k21 + 2 * k22 + 2 * k23 + k24) * h
        y1[i + 1] = y1[i] + (1 / 6) * (k11 + 2 * k12 + 2 * k13 + k14) * h



