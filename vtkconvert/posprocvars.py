import numpy as np
import scipy.io as sio
from scipy.ndimage import gradient
import os

def posprocvars(fstart, fend, p_col, p_row):
    for Flow in range(fstart, fend + 1):
        m = sio.loadmat('mesh.mat')
        
        flowname = f'flow_{Flow:010d}.mat'
        flow = sio.loadmat(flowname)
        U = flow['U']
        U[np.isnan(U)] = 0
        V = flow['V']
        V[np.isnan(V)] = 0
        W = flow['W']
        W[np.isnan(W)] = 0

        Uy, Ux, Uz = np.gradient(U)
        Vy, Vx, Vz = np.gradient(V)
        Wy, Wx, Wz = np.gradient(W)

        Ux = Ux / np.gradient(m['X'].T)
        Uy = Uy / np.gradient(m['Y'])
        Uz = Uz / np.transpose(np.gradient(m['Z']), (0, 2, 1))

        Vx = Vx / np.gradient(m['X'].T)
        Vy = Vy / np.gradient(m['Y'])
        Vz = Vz / np.transpose(np.gradient(m['Z']), (0, 2, 1))

        Wx = Wx / np.gradient(m['X'].T)
        Wy = Wy / np.gradient(m['Y'])
        Wz = Wz / np.transpose(np.gradient(m['Z']), (0, 2, 1))

        L2 = np.zeros_like(U)

        for i in range(U.shape[0]):
            for j in range(U.shape[1] - 60):
                if np.any(U[i, j, :] != 0):
                    for k in range(U.shape[2]):
                        L2[i, j, k] = lambda2(np.array([[Ux[i,j,k], Uy[i,j,k], Uz[i,j,k]],
                                                        [Vx[i,j,k], Vy[i,j,k], Vz[i,j,k]],
                                                        [Wx[i,j,k], Wy[i,j,k], Wz[i,j,k]]]))
                        if L2[i, j, k] >= 0:
                            L2[i, j, k] = np.nan

        q = np.zeros_like(U)
        qs = np.zeros_like(U)

        S11 = Ux
        S12 = 0.5 * (Uy + Vx)
        S13 = 0.5 * (Uz + Wx)
        S22 = Vy
        S23 = 0.5 * (Vz + Wy)
        S33 = Wz
        Omga12 = 0.5 * (Uy - Vx)
        Omga13 = 0.5 * (Uz - Wx)
        Omga23 = 0.5 * (Vz - Wy)

        for i in range(60, U.shape[0]):
            for j in range(U.shape[1] - 60):
                if np.any(U[i, j, :] != 0):
                    for k in range(U.shape[2]):
                        q[i, j, k] = qcrit(Omga12[i,j,k], Omga13[i,j,k], Omga23[i,j,k],
                                           S11[i,j,k], S22[i,j,k], S33[i,j,k],
                                           S12[i,j,k], S13[i,j,k], S23[i,j,k])
                        qs[i, j, k] = np.nan
                        if q[i, j, k] <= 0:
                            qs[i, j, k] = q[i, j, k]
                            q[i, j, k] = np.nan

        mq = np.nanmax(q)
        Q = q / mq
            
        Dilatation = Ux + Vy + Wz

        aux3 = sio.loadmat('mesh.mat')
        gama = 1.4
        X = aux3['X']
        Y = aux3['Y']
        Z = aux3['Z']
        wall = aux3['wall']
        flowParameters = aux3['flowParameters']
        flowType = aux3['flowType']
        # P = 0.4 * flow['R'] * flow['E']
        # T = flow['E'] * gama * (gama - 1) * (flowParameters['Ma']**2)
        Vorti, Vortj, Vortk, Vort = vorticity(flow)
        Mach = (1 / np.sqrt(gama * (gama - 1))) * flow['U'] / np.sqrt(flow['E'])

        sio.savemat(flowname, {'L2': L2, 'Q': Q, 'Dilatation': Dilatation, 'Vortk': Vortk, 'Mach': Mach}, appendmat=True)

def vorticity(flow):
    Vorti, Vortj, Vortk = np.gradient(flow['U'], flow['V'], flow['W'])
    Vort = np.sqrt(Vorti**2 + Vortj**2 + Vortk**2)
    return Vorti, Vortj, Vortk, Vort

def lambda2(J):
    S = (J + J.T) / 2
    O = (J - J.T) / 2
    L = np.linalg.eigvals(np.dot(S, S) + np.dot(O, O))
    return np.real(np.sort(L)[1])

def qcrit(Omga12, Omga13, Omga23, S11, S22, S33, S12, S13, S23):
    A = 2*(Omga12**2 + Omga13**2 + Omga23**2) - (S11**2 + S22**2 + S33**2) - 2*(S12**2 + S13**2 + S23**2)
    return 0.5 * A


