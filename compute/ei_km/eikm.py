from km_util import common, math as mt
from params import prm, storage
import numpy as np
from matplotlib import pyplot as plt
import h5py
from pprint import pprint
import random
from datetime import datetime
import sys


fname = storage + '/' + common.get_filename(storage, 'eiku', '.hdf5', '{:03d}')
print(fname)
input()

pprint(prm)
input()

gammaE = prm['gammaE']
gammaI = prm['gammaI']
omegaE_bar = prm['omegaE_bar'] 
omegaI_bar = prm['omegaI_bar']

KEE = prm['KEE']
KEI = prm['KEI']
KIE = prm['KIE']
KII = prm['KII']

N = prm['N']
r = prm['r']

D = prm['D']

dt = prm['dt']
tmax = prm['tmax']
step = prm['step']

use_cpp = prm['use_cpp']

omegaE = np.array([mt.from_lorentzian(gammaE, omegaE_bar) for i in range(N)]) + KEE - KEI
omegaI = np.array([mt.from_lorentzian(gammaI, omegaI_bar) for i in range(N)]) + KIE - KII

STDEV_NOISE_HALF_TIMESTEP = np.sqrt(4*D/dt)


def ord_param(theta1):
    return np.sum(np.exp(1j*theta1))/N 

def theta_dot(thetaE1, thetaI1):
    zE, zI = ord_param(thetaE1), ord_param(thetaI1)
    rE, rI = np.abs(zE), np.abs(zI) 
    psiE, psiI = np.angle(zE), np.angle(zI)

    return omegaE - (1+r)/2 * (KEE * rE * np.cos(psiE-thetaE1) - KEI * rI * np.cos(psiI-thetaE1)),\
           omegaI - (1+r)/2 * (KIE * rE * np.cos(psiE-thetaI1) - KII * rI * np.cos(psiI-thetaI1))


thetaE = np.array([random.uniform(0, 2*np.pi) for i in range(N)])
thetaI = np.array([random.uniform(0, 2*np.pi) for i in range(N)])

t = 0
i, step = 0, prm['step']
chunk_ind = 0

f = h5py.File(fname, 'w')
f.create_dataset('omegaE', dtype='f', data=omegaE.tolist())
f.create_dataset('omegaI', dtype='f', data=omegaI.tolist())

f.create_dataset('thetaE', dtype='f', data=[thetaE.tolist()], chunks=(1, N), maxshape=(None, N))
f.create_dataset('thetaI', dtype='f', data=[thetaI.tolist()], chunks=(1, N), maxshape=(None, N))

f.create_dataset('t', dtype='f', data=[t], chunks=(1,), maxshape=(None,))

for key in prm:
    f.attrs[key] = prm[key]
f.attrs['time'] = datetime.now().strftime('%Y.%m.%d-%H:%m:%S') 


while t < tmax:
    k1, l1 = theta_dot(thetaE, thetaI)
    k2, l2 = theta_dot(thetaE + dt/2 * k1, thetaI + dt/2 * l1)
    k3, l3 = theta_dot(thetaE + dt/2 * k2, thetaI + dt/2 * l2)
    k4, l4 = theta_dot(thetaE + dt * k3, thetaI + dt * l3)

    
    thetaE = thetaE + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    thetaI = thetaI + dt/6 * (l1 + 2*l2 + 2*l3 + l4)
    t += dt
    i += 1

    if i % step == 0:
        chunk_ind += 1
        sys.stdout.write('\r {:.4f}'.format(t) + f'/{tmax}        ')
        f['thetaE'].resize(chunk_ind+1, axis=0)
        f['thetaE'][chunk_ind] = thetaE.tolist()
        f['thetaI'].resize(chunk_ind+1, axis=0)
        f['thetaI'][chunk_ind] = thetaI.tolist()

        f['t'].resize(chunk_ind+1, axis=0)
        f['t'][chunk_ind] = t

f.close()   
