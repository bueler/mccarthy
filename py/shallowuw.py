"""
This module contains functions which compute the surface values
of the horizontal and vertical velocity using the non-sliding
shallow ice approximation (SIA) in the planar case (one horizontal
dimension).  Inputs are the gridded values of the bed and surface
elevation, as piecewise-linear functions, so analytical gradients
are not required; finite differencing determines gradients.  There
is also a function which computes the diffusivity function.

Plots an example case when run at command line.
"""

import numpy as np

# physical parameters
secpera = 31556925.0
g = 9.81  # acceleration of gravity (m s-2)
rhoi = 910.0  # density of ice (kg m-3)
n = 3  # Glen exponent
A = 3.1689e-24  # EISMINT I value of ice softness (Pa-3 s-1)


def u(x, b, s):
    """Compute surface value of horizontal velocity u_j at
    *interior* regular grid locations.  If input vectors have
    length N then output has length N-2."""
    assert np.all(b <= s)  # check admissibility
    dx = x[1] - x[0]
    dsdx = (s[1:] - s[:-1]) / dx
    Hstag = ((s[:-1] - b[:-1]) + (s[1:] - b[1:])) / 2.0
    C_u = (2.0 / (n + 1)) * A * (rhoi * g) ** n
    ult = -C_u * abs(dsdx[:-1]) ** (n - 1) * dsdx[:-1] * Hstag[:-1] ** (n + 1)
    urt = -C_u * abs(dsdx[1:]) ** (n - 1) * dsdx[1:] * Hstag[1:] ** (n + 1)
    return (ult + urt) / 2.0


def diffusiveD(x, b, s):
    """Compute diffusivity D from SIA, at staggered points.  Used in
    computing the vertical velocity, and to determine a stable time step."""
    dx = x[1] - x[0]
    dsdx = (s[1:] - s[:-1]) / dx
    Hstag = ((s[:-1] - b[:-1]) + (s[1:] - b[1:])) / 2.0
    C_D = (2.0 / (n + 2)) * A * (rhoi * g) ** n
    return C_D * Hstag ** (n + 2) * abs(dsdx) ** (n - 1)


def w(x, b, s):
    """Compute the surface value of the vertical velocity w_j at *interior*
    regular grid locations using centered finite differences.  Applies
    second-order FD formulas to evaluate the divergence of the flux.  If
    input vectors have length N then output has length N-2."""
    dx = x[1] - x[0]
    dsdx = (s[1:] - s[:-1]) / dx
    Hstag = ((s[:-1] - b[:-1]) + (s[1:] - b[1:])) / 2.0
    C_Uds = (2.0 / (n + 1)) * A * (rhoi * g) ** n
    tmp = - C_Uds * (Hstag * abs(dsdx)) ** (n + 1)  # staggered U . ds
    D = diffusiveD(x, b, s)
    return (tmp[:-1] + tmp[1:]) / 2.0 + (D[1:] * dsdx[1:] - D[:-1] * dsdx[:-1]) / dx


# parameters for the plot example
L = 20.0e3  # length of domain (m)
H0 = 400.0  # reference ice thickness (m)
alpha = -0.04  # slope of bed
wave_amp = 25.0  # amplitude of surface waves; set to 0 for slab
wave_len = 8.0e3  # wavelength of surface waves
J = 80  # number of subintervals (in x direction)
outdir = "output/"  # directory name for image output


def b(x):
    "Compute b(x) for the plot example, a slope."
    return alpha * (x - L)


def s(x):
    "Compute s(x) for the plot example, a sine wave on a slope."
    stmp = H0 + wave_amp * np.sin(2.0 * np.pi * x / wave_len)
    return b(x) + np.maximum(stmp, 0.0)


def mkoutdir(dirname):
    try:
        import os

        os.mkdir(dirname)
    except:
        pass


def compute_velocity():
    x = np.linspace(0.0, L, J + 1)
    print(f"  equally-spaced grid:  J = {J} subintervals, dx = {x[1]-x[0]:7.4f} m")
    uu = u(x, b(x), s(x)) * secpera
    ww = w(x, b(x), s(x)) * secpera
    print(f"  range for u:  min = {uu.min():9.4f},  max = {uu.max():9.4f}  (m/a)")
    print(f"  range for w:  min = {ww.min():9.4f},  max = {ww.max():9.4f}  (m/a)")
    if wave_amp == 0.0:
        C_slab = -(2.0 / (n + 1)) * A * (rhoi * g) ** n
        u_slab = C_slab * abs(alpha) ** 2 * alpha * H0 ** (n + 1)
        print(f"  COMPARE: exact slab-on-slope value u={u_slab * secpera:.4f} (m/a)")
    return x, uu, ww


def plot_geometry_velocity(x, u, w):
    import matplotlib.pyplot as plt

    xint = x[1:-1]
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(x / 1.0e3, s(x), ".-", color="C1", ms=4.0, label="s")
    ax1.plot(x / 1.0e3, b(x), "-", color="C0", label="b")
    ax1.legend(loc="lower left")
    ax1.set_xticklabels([])
    ax1.set_ylabel("elevation (m)")
    ax2.plot(xint / 1.0e3, u, "o", ms=6.0, mfc="w", label="u", color="C2")
    ax2.plot(xint / 1.0e3, w, "+", ms=8.0, label="w", color="C3")
    ax2.legend(loc="upper right")
    ax2.set_xticks(np.linspace(0.0, 20.0, 5))
    ax2.set_xticklabels(["0", "5", "10", "15", "20"])
    ax2.set_ylabel("velocity (m/a)")
    ax2.grid(visible=True)
    plt.xlabel("x (km)")
    mkoutdir(outdir)
    print(f"  writing to image file {outdir}uw.pdf")
    plt.savefig(f"{outdir}uw.pdf")


if __name__ == "__main__":
    x, u, w = compute_velocity()
    plot_geometry_velocity(x, u, w)
