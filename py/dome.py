"""
Solves the time-dependent shallow ice approximation model in one
dimension, for s(t,x), under zero surface mass balance a(x)=0.
Surface values of the horizontal u(t,x) and w(t,x) vertical velocity
are computed by shallowuw.py.

Initial surface elevation is either random bumbs or the Halfar (1981)
dome exact solution.

Runs produce .png figures in the output/ directory.

Convenient way to run and view (where "eog" is a .png viewer):
   $ rm -rf output/;  python3 dome.py;  eog output/
"""

import numpy as np
import matplotlib.pyplot as plt
import shallowuw as sia

secpera = 31556925.0

# user parameters
T_a = 0.1  # duration of run in years (a)
L = 20.0e3  # half-length of domain (m)
J = 100  # number of spatial subintervals
halfarinitial = True
c_CFL = 0.005  # coefficient when using CFL to determine time step
maxsteps = 10000  # limit the total number of steps
outdir = "output/"  # directory name for frames

# spatial grid is equally-spaced
x = np.linspace(-L, L, J + 1)  # dx = 2 L / J
dx = x[1] - x[0]

# bed elevation
b = np.zeros(x.shape)

def t0_halfar(x, alpha=1.0/11.0, R0=10.0e3, H0=1400.0):
    # Halfar time-dependent SIA geometry from
    #   P. Halfar (1981), On the dynamics of the ice sheets,
    #   J. Geophys. Res. 86 (C11), 11065--11072
    # For the formula for t0 see equations (3)--(9) in Bueler et al
    # (2005), in 1D and for n=3.
    beta = alpha
    Gamma = 2.0 * sia.A * (sia.rhoi * sia.g)**3.0 / 5.0
    return (7.0 / 4.0)**3.0 * (beta / Gamma) * R0**4.0 / H0**7.0

def s_halfar(t, x, alpha=1.0/11.0, R0=10.0e3, H0=1400.0):
    # Halfar time-dependent SIA geometry from
    #   P. Halfar (1981), On the dynamics of the ice sheets,
    #   J. Geophys. Res. 86 (C11), 11065--11072
    # The solution is evaluated at t = t0.  For the formula for t0 see
    # equations (3)--(9) in Bueler et al (2005), in 1D and for n=3.
    t0 = t0_halfar(x)
    print(f"  evaluating Halfar solution at t = {t / secpera} a, with t0 = {t0 / secpera} a")
    beta = alpha
    pp = 1.0 + 1.0 / 3.0
    rr = 3.0 / (2.0 * 3.0 + 1.0)
    s0 = (t / t0)**beta * R0   # margin position at time t
    xsc = (t / t0)**(-beta) * x[abs(x) < s0] / R0
    s = np.zeros(x.shape)
    s[abs(x) < s0] = H0 * (t / t0)**(-alpha) * (1.0 - abs(xsc)**pp)**rr
    return s

# initial surface elevation
if halfarinitial:
    t0 = t0_halfar(x)
    s = s_halfar(t0, x)
else:
    s = np.zeros(x.shape)
    rng = np.random.default_rng(seed=42)  # repeatable pseudo-random for testing
    bumps = 800.0 * rng.random(x.shape)
    s[abs(x) < L/3] = 800.0 + bumps[abs(x) < L/3]


def skestep(x, s, dt, u, w):
    """Compute a time step of the SKE using the upwind scheme.  The end-point
    values of s are left unmodified.  The input values of u and w are for
    interior points only."""
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

# time-stepping loop
plt.plot(x, s, lw=2.0)
t = 0.0
for k in range(maxsteps):
    #print(".", end="")
    u = sia.u(x, b, s)  # surface horizontal velocity at interior points
    w = sia.w(x, b, s)  # surface vertical velocity at interior points
    dt_CFL = c_CFL * dx / abs(u).max()  # use CFL to determine maximum dt
    dt = min(dt_CFL, T_a * secpera - t)  # do not step past end of [0,T]
    print(f"  t = {t / secpera:.4f} a, step {k}: dt = {dt / secpera:.6f}")
    if dt < 1.0e-6:
        break
    snew = skestep(x, s, dt, u, w)  # take a candidate step
    if np.any(snew < b):
        print(f"    warning: admissibility fails at {sum(snew < b)} locations")
    s = np.maximum(snew, b)  # now force admissibility
    plt.plot(x, s, color="C2")
    t += dt

plt.plot(x, s, lw=2.0, color="C3")
plt.show()
