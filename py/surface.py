"""
Solves time-dependent and steady-state surface kinematical equation
(SKE) models in one dimension, for s(t,x).  The surface mass balance
a(x), horizontal velocity u(x), and vertical velocity w(x) are all
assumed to be given and time-independent.  (This is not realistic!)

Runs produce .png figures in the output/ directory.  Note that a
bed elevation b(x) is shown in these figures, but it does not
enter into the SKE itself, except through admissibility.

Convenient way to run and view (where "eog" is a .png viewer):
   $ rm -rf output/;  python3 surface.py;  eog output/

Explore using code modifications (for example):
   u(x)      -->  2.5 * u(x)  instability
   N = 50    -->  N = 20      instability (dt=1a --> dt=2.5a)
   dx = 1000 -->  dx = 400    instability (J=20 --> J=50)
   a(x)      -->  3 * a(x)    violate s >= b
"""

import numpy as np
import matplotlib.pyplot as plt

secpera = 31556925.0

# major parameters
L = 20.0e3  # length of domain (m)
J = 20  # number of subintervals
T_a = 50  # duration of run (a)
N = 50  # number of time steps
outdir = "output/"  # directory name for frames


def a(x):
    """surface mass balance is decreasing and linear, with
    ELA at x = L/5"""
    return 1.0e-7 - 5.0e-7 * x / L


def b(x):
    """bed elevation is polynomial"""
    z, c2, c3, c4 = L / 2, 4.167e-06, -1.667e-10, -8.333e-15
    return (x - z) ** 2 * (c2 + c3 * (x - z) + c4 * (x - z) ** 2)


def u(x):
    """horizontal velocity is constant"""
    return (500.0 + 0.0 * x) / secpera


def w(x):
    """vertical velocity is zero"""
    return 0.0 * x


def explicitupwindstep(x, s, dt):
    """compute a time step of the SKE, assuming u(x)>0 always,
    by upwind scheme; left-end value snew[0] = s[0] is unmodified"""
    dx = x[1] - x[0]
    snew = s.copy()
    snew[1:] -= dt * u(x[1:]) * (s[1:] - s[:-1]) / dx
    snew[1:] += dt * (a(x[1:]) + w(x[1:]))
    return snew


def mkoutdir(dirname):
    """convenience function for creating an output directory"""
    import os

    try:
        os.mkdir(dirname)
    except FileExistsError:
        pass


def run_evolution():
    # space-time grid info
    x = np.linspace(0.0, L, J + 1)  # dx = (L - 0) / J
    dx = x[1] - x[0]
    dt = (T_a * secpera) / N
    print(f"  stability constant = {u(x).max() * dt / dx:8.4f}")

    # initial surface elevation
    s = b(0.0) - 200.0 * x / L - 50.0 * np.sin(10.0 * x / L)

    # evolution loop: plot frames to .png before each time step
    xplot = np.linspace(0.0, L, 1001)
    plt.figure(1, figsize=(12, 5))
    plt.plot(xplot, b(xplot))
    plt.axis([0, L, 0, 1.5 * b(0.0)])
    plt.xlabel("x (m)")
    plt.ylabel("elevation (m)")
    (sh,) = plt.plot(x, s, "-o")
    th = plt.text(1000.0, 700.0, f"t = {0.0:6.2f} (a)")
    mkoutdir(outdir)
    print(f"  writing {N} steps to image files {outdir}frame???.png")
    for k in range(N + 1):
        plt.savefig(f"{outdir}frame{k:03d}.png")
        s = explicitupwindstep(x, s, dt)
        sh.set_ydata(s)
        th.set_text(f"t = {(k+1) * dt / secpera:6.2f} (a)")


def run_steady():
    # spatial grid
    x = np.linspace(0.0, L, J + 1)
    dx = x[1] - x[0]
    # solve for steady s(x)
    s = x.copy()
    s[0] = b(0.0)
    for j in range(1, J + 1):
        s[j] = s[j - 1] + (dx / u(x[j])) * (a(x[j]) + w(x[j]))
    # save plot to .png
    xplot = np.linspace(0.0, L, 1001)
    plt.figure(2, figsize=(12, 5))
    plt.plot(xplot, b(xplot))
    plt.axis([0, L, 0, 1.5 * b(0.0)])
    plt.xlabel("x (m)")
    plt.ylabel("elevation (m)")
    plt.plot(x, s, "-o")
    plt.text(1000.0, 700.0, "steady state")
    mkoutdir(outdir)
    print(f"  writing to image file {outdir}steady.png")
    plt.savefig(f"{outdir}steady.png")


print("evolution mode ...")
run_evolution()
print("steady state mode ...")
run_steady()
