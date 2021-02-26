import h5py
from km_util import common, hdfutil
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

gamma_glob = 0.4
N_glob = 10000
Kc = 2*gamma_glob

def mask(f: h5py.File):
    return f.attrs['gamma'] == gamma_glob and \
           f.attrs['omega_bar'] == 0 and \
           f.attrs['N'] == N_glob


files = common.get_hdf('data', mask)

K_pl = []
r_pl = []
stdev_pl = []


for f in files:
    attrs = hdfutil.get_attrs(f)

    t = np.array(f['t'])
    theta = np.array(f['theta'])
    N = attrs['N']
    K = attrs['K']
    gamma = attrs['gamma']
    omega_bar = attrs['omega_bar']
    z = np.sum(np.exp(1j*theta), axis=1) / N
    r = np.abs(z) 

    r_ave, r_stdev = hdfutil.get_r(f, 50)

    K_pl.append(K)
    r_pl.append(r_ave)
    stdev_pl.append(r_stdev)

    print(f'K={K}, gamma={gamma}, N={N}, omega_bar={omega_bar}')

plt.scatter(K_pl, r_pl, s=10, marker='s', facecolors='none', edgecolors='r')
plt.errorbar(K_pl, r_pl, yerr=stdev_pl, linestyle="None")

KK = np.linspace(0, max(K_pl), 1000)
plt.plot(KK, np.where(KK > Kc, np.sqrt(1-Kc/KK), 0))
plt.title(f'$\\gamma = {gamma_glob}, N = {N}$')
plt.ylabel('$r$')
plt.xlabel('$K$')
plt.show()
