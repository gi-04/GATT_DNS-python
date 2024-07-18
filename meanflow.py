import os
import scipy.io as sio
import numpy as np

Tin = np.arange(11700, 16601)
pastaIn = 'rugosidade6_5_BC2'

Tout = Tin[-1] + 1
pastaOut = pastaIn

U = V = W = R = E = 0
for t in Tin:
    tstr = os.path.join(pastaIn, f'flow_{t:010d}.mat')
    print(t)
    current = sio.loadmat(tstr)
    if t == Tin[0]:
        U = current['U']
        V = current['V']
        W = current['W']
        R = current['R']
        E = current['E']
    else:
        U += current['U']
        V += current['V']
        W += current['W']
        R += current['R']
        E += current['E']

U /= len(Tin)
V /= len(Tin)
W /= len(Tin)
R /= len(Tin)
E /= len(Tin)
t = current['t']

tstr = os.path.join(pastaOut, f'flow_{Tout:010d}.mat')
sio.savemat(tstr, {'t': t, 'U': U, 'V': V, 'W': W, 'R': R, 'E': E})

    
