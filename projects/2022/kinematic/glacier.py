#!/usr/bin/env python3
# (C) 2022 Ed Bueler

# Generates synthetic, time-dependent glacier geometry, horizontal surface
# motion, and lumped climatic mass balance.  For documentation of formulas
# see synglac.pdf.  For usage instructions see README.md.

import numpy as np

# t in seconds, x in meters, s(t,x) = H(t,x) in meters, u(t,x) in m/s
secpera = 31556926.0
n = 3.0  # Glen exponent
rho = 910.0  # ice density; kg m-3
g = 9.81  # Earth gravity; m s-2
A = 1.0e-16 / secpera  # EISMINT I value
T = 2000.0 * secpera  # 2000 a
Hc_0 = 3000.0  # 3000 m
L_0 = 400.0 * 1.0e3  # 400 km

# derived powers
q = 1.0 + 1.0 / n
r = n / (2.0 * n + 2.0)

# tolerance
icytol = 1.0  # 1.0 m in horizontal

# time-dependent center surface elevation (thickness), and its derivative
def Hc(t):
    return Hc_0 * (1.0 - 0.5 * np.sin(np.pi * t / T))


def dHcdt(t):
    return -0.5 * Hc_0 * np.cos(np.pi * t / T) * (np.pi / T)


# time-dependent glacier length, from center to margin, and its derivative
def L(t):
    return L_0 * (1.0 - 0.75 * np.sin(np.pi * t / T))


def dLdt(t):
    return -0.75 * L_0 * np.cos(np.pi * t / T) * (np.pi / T)


# following functions of (t,x) are vectorized over x but not t

# helper functions
def _geticyy(t, x):
    icy = np.abs(x) < L(t) - icytol
    return icy, (np.abs(x[icy]) / L(t))


def _psi(t, x):
    icy, y = _geticyy(t, x)
    psi = np.zeros(np.shape(x))
    psi[icy] = (n + 1.0) * y - 1.0 + n * (1.0 - y) ** q - n * y**q
    return psi


def _phi(t, x):
    icy, y = _geticyy(t, x)
    phi = np.zeros(np.shape(x))
    phi[icy] = (1.0 - y) ** (1.0 / n) + y ** (1.0 / n) - 1.0
    return phi


def _dpsidt(t, x):
    icy, y = _geticyy(t, x)
    dpsi = np.zeros(np.shape(x))
    dpsi[icy] = (n + 1.0) * (dLdt(t) / L(t)) * y * _phi(t, x)
    return dpsi


def _dpsidx(t, x):
    icy, y = _geticyy(t, x)
    z = np.sign(x[icy]) / L(t)
    dpsi = np.zeros(np.shape(x))
    dpsi[icy] = -(n + 1.0) * z * _phi(t, x)
    return dpsi


# surface elevation s(t,x) from formula (5.50) in van der Veen (2013)
def s(t, x):
    return Hc(t) * (n - 1.0) ** (-r) * _psi(t, x) ** r


# the time rate of change of the surface; dsdt = partial s / partial t
# WARNING: negative power of zero if formula is used where not icy
def dsdt(t, x):
    C = (n - 1.0) ** (-r)
    icy, _ = _geticyy(t, x)
    ds = np.zeros(np.shape(x))
    ds[icy] = C * (
        dHcdt(t) * _psi(t, x[icy]) ** r
        + r * Hc(t) * _psi(t, x[icy]) ** (r - 1.0) * _dpsidt(t, x[icy])
    )
    return ds


# surface slope; dsdx = partial s / partial x
# WARNING: negative power of zero if formula is used where not icy
def dsdx(t, x):
    C = (n - 1.0) ** (-r)
    icy, _ = _geticyy(t, x)
    ds = np.zeros(np.shape(x))
    psi = _psi(t, x[icy])
    assert (psi > 0.0).all(), f"psi supposed to be positive but\n  psi = {psi}"
    ds[icy] = r * Hc(t) * C * psi ** (r - 1.0) * _dpsidx(t, x[icy])
    return ds


# horizontal velocity at surface u(t,x) from SIA formula (see notes)
def us(t, x):
    gam = 2.0 * A * (rho * g) ** n / (n + 1.0)
    return -gam * s(t, x) ** (n + 1.0) * dsdx(t, x) ** n  # assumes n odd integer


# the vertical surface rate atilde(t,x) is the actual SMB a(t,x) lumped
# with the vertical velocity at the surface w|_s(t,x), extended by Inf
def atilde(t, x):
    icy, _ = _geticyy(t, x)
    aa = np.zeros(np.shape(x))
    aa[icy] = dsdt(t, x[icy]) + us(t, x[icy]) * dsdx(t, x[icy])
    # np.NaN would be better missing value but matplotlib mis-handles it
    # (e.g. pcolormesh())
    aa[np.logical_not(icy)] = np.Inf
    return aa


# abalance(t,x) below would *be* the SMB if s(t,x)
# were held at time t; see formula (5.51) in van der Veen (2013)
# [in fact s(t,x) changes in time; atilde(t,x) is the actual SMB,
# but lumped with w|_s]
def abalance(t, x):
    icy = np.abs(x) < L(t)
    Gam = 2.0 * A * (rho * g) ** n / (n + 2.0)
    C = Hc(t) ** (2.0 * n + 2.0) * Gam * (2.0 * L(t) * (1.0 - 1.0 / n)) ** (-n)
    y = np.abs(x[icy]) / L(t)
    phi = y ** (1.0 / n) + (1.0 - y) ** (1.0 / n) - 1.0
    ell = (1.0 - n) / n
    omega = y**ell - (1.0 - y) ** ell
    abal = -(C / L(t)) * np.ones(np.shape(x))
    abal[icy] = (C / L(t)) * phi ** (n - 1.0) * omega
    return abal


# if glacier.py is executed then it produces plots to illustrate
# all important fields
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    print(
        "plotting glacier at t=0, %.0f, %.0f years ..."
        % (0.25 * T / secpera, 0.5 * T / secpera)
    )
    xx = np.linspace(-1.1 * L_0, 1.1 * L_0, 801)

    # figure 1: s at 0,T/4,T/2,3T/4,T
    plt.figure()
    tt = [0.0, 0.25 * T, 0.5 * T, 0.75 * T, T]
    for j in range(5):
        plt.subplot(5, 1, j + 1)
        plt.plot(xx / 1.0e3, s(tt[j], xx), label="%.0f years" % (tt[j] / secpera))
        plt.legend(loc="upper right")
        plt.gca().set_ylim([-200.0, Hc_0 + 200.0])
        plt.gca().set_yticks([0.0, 1500.0, 3000.0])
        if j < 4:
            plt.tick_params(
                axis="x",  # changes apply to the x-axis
                which="both",  # both major and minor ticks are affected
                bottom=False,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                labelbottom=False,
            )  # labels along the bottom edge are off
        else:
            plt.xlabel("x  (km)")
            plt.ylabel("elevation s  (m)")
    # plt.savefig('surfacesnaps.png',bbox_inches='tight')

    # figure 2: s, s_t, s_x in one figure
    plt.figure()
    plt.subplot(3, 1, 1)
    for tt in [0.0, 0.25 * T, 0.5 * T]:
        plt.plot(xx / 1.0e3, s(tt, xx), label="%.0f years" % (tt / secpera))
    plt.legend(loc="upper right")
    plt.ylabel("elevation s  (m)")
    plt.subplot(3, 1, 2)
    for tt in [0.0, 0.25 * T, 0.5 * T]:
        plt.plot(
            xx / 1.0e3,
            dsdt(tt, xx) * secpera,
            ".",
            ms=2.0,
            label="%.0f years" % (tt / secpera),
        )
    plt.ylabel(r"surface rate $\frac{\partial s}{\partial t}$  (m/a)")
    plt.gca().set_ylim([-20.0, 2.0])  # compare range for atilde
    plt.subplot(3, 1, 3)
    for tt in [0.0, 0.25 * T, 0.5 * T]:
        plt.plot(xx / 1.0e3, dsdx(tt, xx), ".", ms=2.0)
    plt.ylabel(r"surface slope $\frac{\partial s}{\partial x}$")
    plt.xlabel("x  (km)")

    # figure 3: us and us*s_x in one figure
    plt.figure()
    plt.subplot(2, 1, 1)
    for tt in [0.0, 0.25 * T, 0.5 * T]:
        plt.plot(xx / 1.0e3, us(tt, xx) * secpera, label="%.0f years" % (tt / secpera))
    plt.legend(loc="lower right")
    plt.ylabel(r"surface velocity $u|_s$  (m/a)")
    plt.subplot(2, 1, 2)
    for tt in [0.0, 0.25 * T, 0.5 * T]:
        plt.plot(xx / 1.0e3, us(tt, xx) * dsdx(tt, xx) * secpera, ".", ms=2.0)
    plt.ylabel(r"product $u|_s \frac{\partial s}{\partial x}$  (m/a)")
    plt.xlabel("x  (km)")

    # figure 4: atilde, abalance
    plt.figure()
    plt.subplot(2, 1, 1)
    for tt in [0.0, 0.25 * T, 0.5 * T]:
        plt.plot(
            xx / 1.0e3,
            atilde(tt, xx) * secpera,
            ".",
            ms=2.0,
            label="%.0f years" % (tt / secpera),
        )
    plt.legend(loc="lower right")
    plt.ylabel(r"$\tilde a$  (m/a)")
    plt.gca().set_ylim([-20.0, 2.0])  # compare range for dsdt
    plt.subplot(2, 1, 2)
    for tt in [0.0, 0.25 * T, 0.5 * T]:
        plt.plot(xx / 1.0e3, abalance(tt, xx) * secpera)
    plt.ylabel("balance a  (m/a)")
    plt.xlabel("x  (km)")

    plt.show()
