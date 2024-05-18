'''
FIXME blurb
Convenient way to run:
   $ rm -rf output/; python3 movesurface.py; eog output/
Exploration:
   0.4 * N etc     instability (N = 50 is threshold)
   2 * u(x) etc    instability
   2 * a(x) etc    violate s >= b
'''

import numpy as np
import matplotlib.pyplot as plt

# major parameters
L = 20.0e3             # length of domain (m)
J = 20                 # number of subintervals
T_a = 100              # duration of run (a)
N = 100                # number of time steps
dname = 'output/'      # directory name for frames

def a(x):
    '''surface mass balance is decreasing and linear'''
    return 1.0e-7 - 5.0e-7 * x / L

def b(x):
    '''bed elevation is polynomial'''
    z, c2, c3, c4 = L / 2, 4.167e-06, -1.667e-10, -8.333e-15
    return (x - z) ** 2 * (c2 + c3 * (x - z) + c4 * (x - z) ** 2)

def u(x):
    '''horizontal velocity is constant'''
    return (500.0 + 0.0 * x) / secpera

def w(x):
    '''vertical velocity is zero'''
    return 0.0 * x

# space-time grid info
x = np.linspace(0.0, L, J + 1)
dx = x[1] - x[0]
secpera = 31556925.0
dt = (T_a * secpera) / N
print(f'stability constant = {u(x).max() * dt / dx:8.4f}')

# initial surface elevation
s = b(0.0) - 200.0 * x / L - 50.0 * np.sin(10.0 * x / L)

def step(s):
    '''compute a time step; note snew[0] = s[0]'''
    snew = s.copy()
    snew[1:] -= dt * u(x[:-1]) * (s[1:] - s[:-1]) / dx
    snew[1:] += dt * (a(x[:-1]) + w(x[:-1]))
    return snew

# evolution loop: plot frames to .png before each time step
xplot = np.linspace(0.0, L, 1001)
plt.plot(xplot, b(xplot))
plt.axis([0, L, 0, 1.5 * b(0.0)])
plt.xlabel('x (m)')
plt.ylabel('elevation (m)')
sh, = plt.plot(x, s, "-o")
th = plt.text(1000.0, 700.0, f't = {0.0:6.2f} (a)', name='DejaVu Sans Mono')
import os
try:
    os.mkdir(dname)
except FileExistsError:
    pass
for k in range(N+1):
    plt.savefig(f'{dname}frame{k:03d}.png')
    s = step(s) # take a time step
    sh.set_ydata(s)
    th.set_text(f't = {(k+1) * dt / secpera:6.2f} (a)')
