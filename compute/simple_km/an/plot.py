import h5py
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from km_util import common

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

filename = 'ku010.hdf5'

f = h5py.File(f'data/{filename}', 'r')
omega = np.array(f['omega'])
theta = np.array(f['theta'])
t = np.array(f['t'])

gamma = f.attrs['gamma']
N = f.attrs['N']
K = f.attrs['K']

z = np.sum(np.exp(1j*theta), axis=1) / N
theta_ave = np.sum(theta, axis=1) / N
print(f'N = {f.attrs["N"]}')
print(f'K = {f.attrs["K"]}')
print(f'gam = {f.attrs["gamma"]}')

plt.scatter(t, np.abs(z), s=2)
plt.xlabel('$t$', fontsize=18)
plt.ylabel('$r$', fontsize=18)
plt.title(f'$\\gamma = {gamma}, K = {K}, N = {N}$', fontsize=20)
plt.show()

plt.scatter(t, np.sin(np.angle(z)), s=2)
plt.show()

plt.scatter(t, theta_ave, s=2)

plt.show()

f.close()
