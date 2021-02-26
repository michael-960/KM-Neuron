import h5py
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from km_util import common

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

filename = 'eiku000.hdf5'

f = h5py.File(f'data/{filename}', 'r')
omegaE = np.array(f['omegaE'])
omegaI = np.array(f['omegaI'])
thetaE = np.array(f['thetaE'])
thetaI = np.array(f['thetaI'])

t = np.array(f['t'])

omegaE_bar = f.attrs['omegaE_bar']
omegaI_bar = f.attrs['omegaI_bar']
gammaE = f.attrs['gammaE']
gammaI = f.attrs['gammaI']

N = f.attrs['N']
r = f.attrs['r']
rr = min(0.98, r)

KEE = f.attrs['KEE']
KEI = f.attrs['KEI']
KIE = f.attrs['KIE']
KII = f.attrs['KII']

zE = np.sum(np.exp(1j*thetaE), axis=1) / N
zI = np.sum(np.exp(1j*thetaI), axis=1) / N
rE = np.abs(zE)
rI = np.abs(zI)
psiE = np.angle(zE)
psiI = np.angle(zI)

hE = (1-rE**2)/(1+rE**2-2*rE*np.cos(psiE))
hI = (1-rI**2)/(1+rI**2-2*rI*np.cos(psiI))

thetaE_ave = np.sum(thetaE, axis=1) / N
print(f'N = {N}')
print(f'r = {r}')
print(f'gamE = {gammaE}')
print(f'gamI = {gammaI}')


plt.scatter(t, rE, s=0.5, color='red', label='$R_E$')
plt.scatter(t, rI, s=0.5, color='blue', label='$R_I$')

plt.xlabel('$t$', fontsize=18)
plt.ylabel('$R_E, R_I$', fontsize=18)
plt.title(f'$\\bar\\omega_E = {omegaE_bar}, \\bar\\omega_I = {omegaI_bar}, R = {r}, N = {N}$', fontsize=20)
plt.ylim(0, 0.8)
plt.legend()
plt.show()

plt.plot(t, hE, color='red', label='$h_E$')
plt.plot(t, hI, color='blue', label='$h_I$')
plt.ylim(0, 10)
plt.xlabel('$t$', fontsize=18)
plt.ylabel('$h_E, h_I$', fontsize=18)
plt.title(f'$\\bar\\omega_E = {omegaE_bar}, \\bar\\omega_I = {omegaI_bar}, R = {r}, N = {N}$', fontsize=20)
plt.legend()
plt.show()

