"""
Solves the time-dependent shallow ice approximation model in one dimension,
for s(t,x), under zero surface mass balance a(x)=0.  The dynamics model is
nonsliding, isothermal SIA with n=3, on a flat bed.  Initial surface
elevation is either random bumps (default) or the Halfar (1981) dome exact
solution.  Surface values of the horizontal u(t,x) and vertical w(t,x)
velocity are computed by shallowuw.py.  Applies adaptive time-stepping
based on diffusive and CFL conditions.

Runs produce two .png figures in the output/ directory.  Here is a
convenient way to run and view (where "eog" is an image viewer):
   $ python3 dome.py;  eog output/
"""

import numpy as np
import matplotlib.pyplot as plt
import shallowuw as sia

# choose initial condition
halfarinitial = False

# user parameters
T_a = 40.0  # duration of run in years (a)
L = 20.0e3  # half-length of domain (m)
J = 100  # number of spatial subintervals
maxsteps = 1000000  # limit the total number of steps
outdir = "output/"  # directory name for images

secpera = 31556925.0

# spatial grid is equally-spaced
x = np.linspace(-L, L, J + 1)  # dx = 2 L / J
dx = x[1] - x[0]

# bed elevation
b = np.zeros(x.shape)


def t0_halfar(x, alpha=1.0 / 11.0, R0=10.0e3, H0=1200.0):
    # Halfar time-dependent SIA geometry from
    #   P. Halfar (1981), On the dynamics of the ice sheets,
    #   J. Geophys. Res. 86 (C11), 11065--11072
    # For the formula for t0 see equations (3)--(9) in Bueler et al
    # (2005), in 1D and for n=3.
    beta = alpha
    Gamma = 2.0 * sia.A * (sia.rhoi * sia.g) ** 3.0 / 5.0
    return (7.0 / 4.0) ** 3.0 * (beta / Gamma) * R0**4.0 / H0**7.0


def s_halfar(x, alpha=1.0 / 11.0, R0=10.0e3, H0=1200.0):
    # Halfar time-dependent SIA geometry from
    #   P. Halfar (1981), On the dynamics of the ice sheets,
    #   J. Geophys. Res. 86 (C11), 11065--11072
    # The solution is evaluated at t = t0.  For the formula for t0 see
    # equations (3)--(9) in Bueler et al (2005), in 1D and for n=3.
    t0 = t0_halfar(x)
    print(f"  [Halfar solution has t0 = {t0 / secpera:.4f} a]")
    beta = alpha
    pp = 1.0 + 1.0 / 3.0
    rr = 3.0 / (2.0 * 3.0 + 1.0)
    t = t0  # change this to different ratios  t = c t0  to see Halfar solution behavior
    s0 = (t / t0) ** beta * R0  # margin position at time t
    xsc = (t / t0) ** (-beta) * x[abs(x) < s0] / R0
    s = np.zeros(x.shape)
    s[abs(x) < s0] = H0 * (t / t0) ** (-alpha) * (1.0 - abs(xsc) ** pp) ** rr
    return s


# initial surface elevation
if halfarinitial:
    s = s_halfar(x)
else:
    s = np.zeros(x.shape)
    # make following repeatable pseudo-random by setting seed (for testing)
    rng = np.random.default_rng()
    bumps = 800.0 * rng.random(x.shape)
    s[abs(x) < L / 3] = 800.0 + bumps[abs(x) < L / 3]


def getdt(t, tf, maxu, maxD):
    """Determine time step using CFL and diffusive stability conditions,
    and final time.  Returns a very conservative time step."""
    c_D = 1.0 / 6.0  # coefficient from Hindmarsh (2001), for n=3
    dt_D = c_D * dx**2 / maxD  # apply diffusive condition
    c_CFL = 0.9  # go close to CFL limit
    dt_CFL = c_CFL * dx / maxu  # apply CFL condition
    dt_end = tf - t  # do not step past end of [0,tf]
    dts = np.array([dt_D, dt_CFL, dt_end])
    return min(dts), np.argmin(dts)


def skestep(x, s, dt, u, w):
    """Compute a time step of the SKE using the upwind scheme in space
    and forward Euler in time.  The end-point values of s are left unmodified.
    The input values of u and w are for interior points only.  The time step
    dt must be set before calling this."""
    dx = x[1] - x[0]
    J = len(x) - 1
    assert len(s) == J + 1 and len(u) == J - 1 and len(w) == J - 1
    snew = s.copy()
    for j in range(J - 1):
        if u[j] > 0.0:
            snew[j + 1] -= dt * u[j] * (s[j + 1] - s[j]) / dx
        else:
            snew[j + 1] -= dt * u[j] * (s[j + 2] - s[j + 1]) / dx
        snew[j + 1] += dt * w[j]
    return snew


# time-stepping loop
plt.figure
plt.plot(x / 1000.0, s, label="initial: t = 0.00 a")
t = 0.0
tf = T_a * secpera
timedata = []
print("starting time-stepping ...")
for k in range(maxsteps):
    u = sia.u(x, b, s)  # surface horizontal velocity at interior points
    w = sia.w(x, b, s)  # surface vertical velocity at interior points
    maxD = sia.diffusiveD(x, b, s).max()
    dt, reason = getdt(t, tf, abs(u).max(), maxD)
    timedata.append([t, dt, reason])
    # print(f"  step {k}: t = {t / secpera:.4f} a, dt = {dt / secpera:.3e} a")
    if tf - t < 1.0e-6:  # normal stopping criterion
        break
    if dt < 1.0e-6:  # time step less than millionth of second is small
        print(f"HALTING because time step is too small")
        break
    snew = skestep(x, s, dt, u, w)  # take a candidate step
    if np.any(snew < b):
        print(f"[warning: admissibility fails at {sum(snew < b)} locations]")
    s = np.maximum(snew, b)  # now force admissibility
    plt.plot(x / 1000.0, s, lw=0.2, color="C2", alpha=0.2)
    t += dt
print(f"took {k} steps to reach final time {T_a:.4f} a")

plt.plot(x / 1000.0, s, color="C3", label=f"final: t={t / secpera:.2f} a")
plt.xlabel("x (km)")
plt.ylabel("elevation (m)")
plt.legend()
# plt.show()
sia.mkoutdir(outdir)
print(f"writing simulation result to {outdir}dome.png")
plt.savefig(f"{outdir}dome.png")
plt.close()

plt.figure
td = np.array(timedata[:-1])
tdD = td[td[:, 2] == 0, 0:2] / secpera
plt.semilogy(tdD[:, 0], tdD[:, 1], ".", color="C1", label="diffusive")
tdC = td[td[:, 2] == 1, 0:2] / secpera
plt.semilogy(tdC[:, 0], tdC[:, 1], ".", color="C2", label="CFL")
tdE = td[td[:, 2] == 2, 0:2] / secpera
plt.semilogy(tdE[:, 0], tdE[:, 1], ".", color="C4", label="final time")
plt.legend()
plt.grid("on")
plt.xlabel("t (a)")
plt.ylabel("time step (a)")
# plt.show()
print(f"writing time steps graph to {outdir}timesteps.png")
plt.savefig(f"{outdir}timesteps.png")
plt.close()
