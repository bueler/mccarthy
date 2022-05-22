#!/usr/bin/env python3
# (C) 2022 Ed Bueler

import numpy as np

# t in seconds, x in meters, s(t,x) = H(t,x) in meters, u(t,x) in m/s
secpera = 31556926.0
n = 3.0
rho = 910.0
g = 9.81
A = 1.0e-16/secpera    # EISMINT I value
T = 2000.0 * secpera   # 2000 a
H0_0 = 3000.0          # 3000 m
L_0 = 400.0 * 1.0e3    # 400 km

# derived constants
q = 1.0 + 1.0 / n
r = n / (2.0 * n + 2.0)

# time-dependent center surface elevation (thickness), and its derivative
def H0(t):
    return H0_0 * (1.0 - 0.5 * np.sin(np.pi * t / T))

def dH0dt(t):
    return - 0.5 * H0_0 * np.cos(np.pi * t / T) * (np.pi / T)

# time-dependent glacier length, from center to margin, and its derivative
def L(t):
    return L_0 * (1.0 - 0.75 * np.sin(np.pi * t / T))

def dLdt(t):
    return - 0.75 * L_0 * np.cos(np.pi * t / T) * (np.pi / T)

# following functions of (t,x) are vectorized over x but not t

# helper function psi(t,x) and its time derivative
def _psi(t,x):
    icy = (np.abs(x) < L(t))
    y = np.abs(x[icy]) / L(t)
    psi = np.zeros(np.shape(x))
    psi[icy] = (n+1.0) * y - 1.0 + n * (1.0 - y)**q - n * y**q
    return psi

def _dpsidt(t,x):
    icy = (np.abs(x) < L(t))
    y = np.abs(x[icy]) / L(t)
    phi = y**(1.0/n) + (1.0 - y)**(1.0/n) - 1.0
    dpsi = np.zeros(np.shape(x))
    dpsi[icy] = (n+1.0) * y * (dLdt(t) / L(t)) * phi
    return dpsi

# surface elevation s(t,x) from formula (5.50) in van der Veen (2013)
def s(t,x):
    C = H0(t) / (n - 1.0)**r
    return C * _psi(t,x)**r

# surface slope ds/dx(t,x) = partial s / partial t
def dsdx(t,x):
    icy = (np.abs(x) < L(t))
    C = r * H0(t) / (n - 1.0)**r
    gamma = np.zeros(np.shape(x))
    gamma[icy] = _psi(t,x[icy])**(r-1.0)
    omega = np.zeros(np.shape(x))
    y = np.abs(x[icy]) / L(t)
    z = np.sign(x[icy]) / L(t)
    omega[icy] = (n+1.0) * z - (n+1.0) * (1.0 - y)**(1.0/n) * z \
                             - (n+1.0) * y**(1.0/n)  * z
    return C * gamma * omega

# horizontal velocity at surface u(t,x) from SIA formula (see notes)
def us(t,x):
    gam = 2.0 * A * (rho * g)**n / (n+1.0)
    return - gam * s(t,x)**(n+1.0) * dsdx(t,x)**n   # assumes n odd integer

# the surface mass balance a_balance(t,x) would *be* the SMB if s(t,x)
# were held at time t; see formula (5.51) in van der Veen (2013)
# [in fact s(t,x) changes in time; atrue(t,x) is the actual SMB]
def abalance(t,x):
    icy = (np.abs(x) < L(t))
    Gam = 2.0 * A * (rho * g)**n / (n+2.0)
    C = H0(t)**(2.0*n+2.0) * Gam * (2.0 * L(t) * (1.0 - 1.0/n))**(-n)
    y = np.abs(x[icy]) / L(t)
    phi = y**(1.0/n) + (1.0 - y)**(1.0/n) - 1.0
    ell = (1.0-n) / n
    omega = y**ell - (1.0 - y)**ell
    abal = - (C / L(t)) * np.ones(np.shape(x))
    abal[icy] = (C / L(t)) * phi**(n-1.0) * omega
    return abal

# the time rate of change of the surface; \partial s / \partial t
def dsdt(t,x):
    C = (n-1.0)**(-r)
    icy = (np.abs(x) < L(t))
    ds = np.zeros(np.shape(x))
    ds[icy] = C * ( dH0dt(t) * _psi(t,x[icy])**r \
                   + r * H0(t) * _psi(t,x[icy])**(r-1.0) * _dpsidt(t,x[icy]) )
    return ds

# the surface mass balance atrue(t,x) is the actual SMB
# FIXME extend by ablation I think
def atrue(t,x):
    # FIXME
    return np.zeros(np.shape(x))

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    print('plotting glacier at t=0, %.0f, %.0f years ...' \
          % (0.25*T/secpera,0.5*T/secpera))
    xx = np.linspace(-1.1*L_0,1.1*L_0,801)

    # figure 1: s, s_x, u
    plt.figure()
    plt.subplot(3,1,1)
    for tt in [0.0, 0.25*T, 0.5*T]:
        plt.plot(xx/1.0e3,s(tt,xx),
                 label='%.0f years' % (tt/secpera))
    plt.ylabel('elevation s  (m)')
    plt.legend()
    plt.subplot(3,1,2)
    for tt in [0.0, 0.25*T, 0.5*T]:
        plt.plot(xx/1.0e3,dsdx(tt,xx),'.',ms=2.0)
    plt.ylabel('surface slope s_x')
    plt.subplot(3,1,3)
    for tt in [0.0, 0.25*T, 0.5*T]:
        plt.plot(xx/1.0e3,us(tt,xx)*secpera)
    plt.ylabel('surface velocity u  (m/a)')
    plt.xlabel('x  (km)')

    # figure 2: s_t, abalance, atrue
    plt.figure()
    plt.subplot(3,1,1)
    for tt in [0.0, 0.25*T, 0.5*T]:
        plt.plot(xx/1.0e3,dsdt(tt,xx)*secpera,'.',ms=2.0,
                 label='%.0f years' % (tt/secpera))
    plt.ylabel('surface rate s_t  (m/a)')
    plt.legend()
    plt.subplot(3,1,2)
    for tt in [0.0, 0.25*T, 0.5*T]:
        plt.plot(xx/1.0e3,abalance(tt,xx)*secpera)
    plt.ylabel('balance a  (m/a)')
    plt.subplot(3,1,3)
    for tt in [0.0, 0.25*T, 0.5*T]:
        plt.plot(xx/1.0e3,atrue(tt,xx)*secpera)
    plt.ylabel('true a  (m/a)')
    plt.xlabel('x  (km)')

    plt.show()
