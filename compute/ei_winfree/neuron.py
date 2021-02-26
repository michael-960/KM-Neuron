import random
import numpy as np
from matplotlib import pyplot as plt
from pprint import pprint
import sys
import os
from datetime import datetime

import params as prm


N = prm.N
gamma = prm.gamma
r = prm.r
use_cpp = prm.use_cpp
tmpfile_key = f'{datetime.now()}'.replace(' ', '-').replace(':', '-').replace('.', '-')
tmpi = f'tmp/i{tmpfile_key}.tmp'
tmpo = f'tmp/o{tmpfile_key}.tmp'

#%% prepare omega
omega_bar_E = prm.omega_bar_E
omega_bar_I = prm.omega_bar_I
K_EE = prm.K_EE
K_II = prm.K_II
K_EI = prm.K_EI
K_IE = prm.K_IE


def omega_from_lorentzian(gamma, omega_bar):
    y = random.uniform(0, 1)
    omega = gamma * np.tan(np.pi*(y-1/2)) + omega_bar
    return omega

omega_E = []
omega_I = []
for i in range(N):
    omega_E.append(omega_from_lorentzian(gamma, omega_bar_E))
    omega_I.append(omega_from_lorentzian(gamma, omega_bar_I))
    
omega_E = np.array(omega_E)
omega_I = np.array(omega_I)

#%% prepare theta
theta_E = np.array([random.uniform(0, np.pi*2) for i in range(N)])
theta_I = np.array([random.uniform(0, np.pi*2) for i in range(N)])

#%% theta_dot
def P(theta):
    return (1-r)*(1+np.cos(theta))/(1-2*r*np.cos(theta)+r**2)

def Q(theta):
    return 1 - np.cos(theta)

def h(theta_E1, theta_I1):
    return 1/N * np.sum(P(theta_E1)), 1/N * np.sum(P(theta_I1))

def theta_dot(theta_E1, theta_I1):
    h_E, h_I = h(theta_E1, theta_I1)
    return omega_E + Q(theta_E1) * (K_EE * h_E - K_EI * h_I), omega_I + Q(theta_I1) * (K_IE * h_E - K_II * h_I)

#%% simulation
dt = prm.dt
t_max = prm.t_max

t = 0
i = 0
step = 17

tr = []
hEr = []
hIr = []

if use_cpp:
    print(f'using cpp, temp file key = {tmpfile_key}')
    dat = '' #N, K_EE, K_II, K_EI, K_IE, r, theta_E, theta_I, omega_E, omega_I, t_max, dt
    dat = f'{N}\n{K_EE}\n{K_II}\n{K_EI}\n{K_IE}\n{r}\n{t_max}\n{dt}\n\n\n'
    for theta in theta_E:
        dat += f'{theta}\n'
    dat += '\n'
    for theta in theta_I:
        dat += f'{theta}\n'
    dat += '\n'
    for omega in omega_E:
        dat += f'{omega}\n'
    dat += '\n'
    for omega in omega_I:
        dat += f'{omega}\n'
    dat += '\n'
    with open(tmpi, 'w+') as f:
        f.write(dat)

    os.system(f'./neuron {tmpfile_key}')
    
    with open(tmpo, 'r') as f:
        dat_s = f.read()
    dat1 = dat_s.split('/')
    datt = dat1[0].split()
    dathE = dat1[1].split()
    dathI = dat1[2].split()
    for k in range(len(datt)):
        tr.append(float(datt[k]))
        hEr.append(float(dathE[k]))
        hIr.append(float(dathI[k]))
else:
    while t < t_max:
        kE1, kI1 = theta_dot(theta_E, theta_I)
        kE2, kI2 = theta_dot(theta_E + dt/2 * kE1, theta_I + dt/2 * kI1)
        kE3, kI3 = theta_dot(theta_E + dt/2 * kE2, theta_I + dt/2 * kI2)
        kE4, kI4 = theta_dot(theta_E + dt * kE3, theta_I + dt * kI3)
        
        theta_E = theta_E + dt/6 * (kE1 + 2*kE2 + 2*kE3 + kE4)
        theta_I = theta_I + dt/6 * (kI1 + 2*kI2 + 2*kI3 + kI4)
        t += dt
        i += 1

        if i % step == 0:
            sys.stdout.write(f'\r {t}     ')
            tr.append(t)
            h_E, h_I = h(theta_E, theta_I)
            hEr.append(h_E)
            hIr.append(h_I)

print()
tr = np.array(tr)
hEr = np.array(hEr)
hIr = np.array(hIr)

plt.scatter(tr, hEr, s=1)
plt.scatter(tr, hIr, s=1)
plt.xlabel('t')
plt.ylabel('$h_E, h_I$')
plt.show()

