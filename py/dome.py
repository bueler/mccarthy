"""
Solves the time-dependent shallow ice approximation model in one dimension,
for s(t,x), under zero surface mass balance a(x)=0.  The dynamics model is
nonsliding, isothermal SIA with n=3, on a flat bed.  The SIA diffusivity is
computed by a function in shallowuw.py.  Applies adaptive time-stepping
based on the diffusive stability condition.

The initial surface elevation is either random bumps (the default) or the
Halfar (1981) dome exact solution.  In the later case we report the numerical
error in surface elevation.

Runs produce two .png figures in the output/ directory.  Here is a
convenient way to run and view (where "eog" is an image viewer):
   $ python3 dome.py;  eog output/
"""

import numpy as np
import matplotlib.pyplot as plt
import shallow as sia

# choose initial condition
halfarinitial = False

# user parameters
T_a = 40.0  # duration of run in years (a)
L = 20.0e3  # half-length of domain (m)
J = 80  # number of spatial subintervals
maxsteps = 1000000  # limit the total number of steps
outdir = "output/"  # directory name for images

secpera = 31556925.0

# spatial grid is equally-spaced
x = np.linspace(-L, L, J + 1)  # dx = 2 L / J
dx = x[1] - x[0]

# bed elevation
b = np.zeros(x.shape)


def t0_halfar(alpha=1.0 / 11.0, R0=10.0e3, H0=1200.0):
    # For the formula for t0 see equations (3)--(9) in
    #   E. Bueler, C. S. Lingle, J. A. Kallen-Brown, D. N. Covey,
    #   & L. N. Bowman (2005). Exact solutions and verification of
    #   numerical models for isothermal ice sheets. J. Glaciol. 51
    #   (173), 291--306
    beta = alpha
    Gamma = 2.0 * sia.A * (sia.rhoi * sia.g) ** 3.0 / 5.0
    return (7.0 / 4.0) ** 3.0 * (beta / Gamma) * R0**4.0 / H0**7.0


def s_halfar(t, x, alpha=1.0 / 11.0, R0=10.0e3, H0=1200.0):
    # Halfar time-dependent SIA geometry (surface elevation) from
    #   P. Halfar (1981). On the dynamics of the ice sheets,
    #   J. Geophys. Res. 86 (C11), 11065--11072
    # This exact solution is evaluated at an arbitrary time t.
    t0 = t0_halfar()
    beta = alpha
    pp = 1.0 + 1.0 / 3.0
    rr = 3.0 / (2.0 * 3.0 + 1.0)
    s0 = (t / t0) ** beta * R0  # margin position at time t
    xsc = (t / t0) ** (-beta) * x[abs(x) < s0] / R0
    s = np.zeros(x.shape)
    s[abs(x) < s0] = H0 * (t / t0) ** (-alpha) * (1.0 - abs(xsc) ** pp) ** rr
    return s


# initial surface elevation
if halfarinitial:
    t0h = t0_halfar()
    s = s_halfar(t0h, x)
    print(f"  [initial state is Halfar solution at t0 = {t0h / secpera:.4f} a]")
else:
    # put random bumps in middle third of domain
    s = np.zeros(x.shape)
    rng = np.random.default_rng()  # set seed for repeatable pseudo-random
    bumps = 800.0 * rng.random(x.shape)
    s[abs(x) < L / 3] = 800.0 + bumps[abs(x) < L / 3]


def getdt(t, tf, maxD):
    """Determine stable time step using diffusive stability condition,
    and final time."""
    c_D = 1.0 / 6.0  # coefficient from Hindmarsh (2001), for n=3
    dt_D = c_D * dx**2 / maxD  # apply diffusive condition
    return min([dt_D, tf - t])  # do not step past end of [0,tf]


def geometrystep(x, s, D, dt):
    """Compute a time step of the SKE, in its planar a=0 diffusion form
        s_t = (D s_x)_x,
    using forward Euler in time and the centered-difference in space.
    The end-point values of s are left unmodified.  The time step
    dt must be set before calling this."""
    assert len(x) == len(s)
    assert len(D) == len(x) - 1
    dx = x[1] - x[0]
    dsdx = (s[1:] - s[:-1]) / dx
    snew = s.copy()
    snew[1:-1] += dt * (D[1:] * dsdx[1:] - D[:-1] * dsdx[:-1]) / dx
    return snew


# time-stepping loop
plt.figure
plt.plot(x / 1000.0, s, label="initial: t = 0.00 a")
t = 0.0
tf = T_a * secpera
timedata = []
print(f"grid: J={J} intervals, spacing dx={dx:.1f} m")
print(f"starting time-stepping ...")
for k in range(maxsteps):
    D = sia.diffusiveD(x, b, s)
    dt = getdt(t, tf, D.max())
    timedata.append([t, dt])
    if tf - t < 1.0e-6:  # normal stopping criterion
        break
    if dt < 1.0e-6:  # time step less than millionth of second is small
        print(f"HALTING EARLY because time step is too small")
        break
    s = geometrystep(x, s, D, dt)  # take a candidate step
    #if np.any(s < b):
    #    print(f"[admissibility fails at t={t / secpera}, {sum(snew < b)} locations]")
    s = np.maximum(s, b)  # force admissibility
    plt.plot(x / 1000.0, s, lw=0.2, color="C2", alpha=0.2)
    t += dt
print(f"took {k} steps to reach final time {T_a:.4f} a")

plt.plot(x / 1000.0, s, color="C3", label=f"final: t={t / secpera:.2f} a")
if halfarinitial:
    sexact = s_halfar(t0h + tf, x)
    serr = abs(s - sexact)
    serrav = sum(serr) / len(serr[serr > 0.0])
    print(f"  [errors at final time: |s - sexact|_max = {serr.max():.2f} m, |s - sexact|_av = {serrav:.2f} m")
    plt.plot(x / 1000.0, sexact, '.', color="C4", label=f"final exact")

plt.xlabel("x (km)")
plt.ylabel("elevation (m)")
plt.legend()
# plt.show()
sia.mkoutdir(outdir)
print(f"writing results to {outdir}dome.png and {outdir}timesteps.png")
plt.savefig(f"{outdir}dome.png")
plt.close()

plt.figure
td = np.array(timedata[:-1]) / secpera
plt.semilogy(td[:, 0], td[:, 1], ".")
plt.grid("on")
plt.axis([-0.01*T_a, 1.01*T_a, td[:, 1].min(), 1.0])
plt.xlabel("t (a)")
plt.ylabel("time step (a)")
# plt.show()
plt.savefig(f"{outdir}timesteps.png")
plt.close()
