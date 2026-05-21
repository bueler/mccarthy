"""
Solves time-dependent shallow ice approximation model in one
dimension, for s(t,x), under zero surface mass balance a(x)=0.
Surface values of horizontal velocity u(t,x) and vertical velocity
w(t,x) are computed by shallowuw.py.

Runs produce .png figures in the output/ directory.

Convenient way to run and view (where "eog" is a .png viewer):
   $ rm -rf output/;  python3 dome1d.py;  eog output/
"""

import numpy as np
import matplotlib.pyplot as plt
import shallowuw as sia

secpera = 31556925.0

L = 20.0e3  # half-length of domain (m)
J = 40  # number of subintervals
T_a = 0.001  # duration of run (a)
N = 200  # number of time steps
outdir = "output/"  # directory name for frames

x = np.linspace(-L, L, J + 1)  # dx = 2 L / J
dx = x[1] - x[0]
dt = (T_a * secpera) / N

s = np.zeros(x.shape)
rng = np.random.default_rng()
bumps = 700.0 * rng.random(x.shape)
s[abs(x) < L/3] = 700.0 + bumps[abs(x) < L/3]

b = np.zeros(x.shape)

#plt.plot(x, s)
#plt.show()
    
def skestep(x, s, dt, u, w):
    """Compute a time step of the SKE using upwind scheme.  The end-point
    values of s are left unmodified."""
    dx = x[1] - x[0]
    J = len(x) - 1
    assert len(s) == J + 1 and len(u) == J - 1 and len(w) == J - 1
    snew = s.copy()
    for j in range(J - 1):
        if u[j] > 0.0:
            snew[j+1] -= dt * u[j] * (s[j+1] - s[j]) / dx
        else:
            snew[j+1] -= dt * u[j] * (s[j+2] - s[j+1]) / dx
        snew[j+1] += dt * w[j]
    return snew

plt.plot(x, s, lw=2.0)
for k in range(N):
    print(".", end="")
    u = sia.u(x, b, s)  # horizontal velocity at interior points
    w = sia.w(x, b, s)  # vertical velocity at interior points
    snew = skestep(x, s, dt, u, w)
    if np.any(snew < b):
        print(f"\n  [admissibility fails at {sum(snew < b)} locations]")
    s = np.maximum(snew, b)
    plt.plot(x, s, color="C2")

print("\n")
plt.show()
