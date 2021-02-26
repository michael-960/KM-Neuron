from km_util import math as mt
from km_util import common
import numpy as np
from pprint import pprint
import random
import sys
import scipy, numpy
import h5py
import ujson
from params import prm, storage
from datetime import datetime

fname = storage + '/' + common.get_filename(storage, 'ku', '.hdf5', '{:03d}')
print(fname)
input()

pprint(prm)
input()

gamma = prm['gamma']
omega_bar = prm['omega_bar']
K = prm['K']
N = prm['N']
dt = prm['dt']
tmax = prm['tmax']


omega = np.array([mt.from_lorentzian(gamma, omega_bar) for i in range(N)])

def ord_param(theta1):
    S = np.sum(np.exp(1j*theta1)) / N
    return np.abs(S), np.angle(S)

def theta_dot(theta1):
    r, psi = ord_param(theta1)
    return omega + K*r*np.sin(psi-theta1)

theta = np.array([random.uniform(0, 2*np.pi) for i in range(N)])

t = 0
i, step = 0, prm['step']
chunk_ind = 0

f = h5py.File(fname, 'w')
f.create_dataset('omega', dtype='f', data=omega.tolist())
f.create_dataset('theta', dtype='f', data=[theta.tolist()], chunks=(1, N), maxshape=(None, N))
f.create_dataset('t', dtype='f', data=[t], chunks=(1,), maxshape=(None,))

for key in prm:
    f.attrs[key] = prm[key]
f.attrs['time'] = datetime.now().strftime('%Y.%m.%d-%H:%m:%S') 

while t < tmax:
    k1 = theta_dot(theta)
    k2 = theta_dot(theta + dt/2 * k1)
    k3 = theta_dot(theta + dt/2 * k2)
    k4 = theta_dot(theta + dt * k3)
    
    theta = theta + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    t += dt
    i += 1

    if i % step == 0:
        chunk_ind += 1
        sys.stdout.write('\r {:.4f}'.format(t) + f'/{tmax}        ')
        f['theta'].resize(chunk_ind+1, axis=0)
        f['theta'][chunk_ind] = theta.tolist()
        f['t'].resize(chunk_ind+1, axis=0)
        f['t'][chunk_ind] = t

f.close()
