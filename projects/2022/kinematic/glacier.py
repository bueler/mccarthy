#!/usr/bin/env python3
# (C) 2022 Ed Bueler

import numpy as np

# t in seconds, x in meters, s(t,x) = H(t,x) in meters, u(t,x) in m/s
secpera = 31556926.0
n = 3.0
rho = 910.0
g = 9.81
A = 1.0e-16/secpera    # EISMINT I value
_T = 1000.0 * secpera  # 1000 a
_H0 = 1000.0
_L = 400.0 * 1.0e3     # 400 km

# derived constants
_q = 1.0 + 1.0 / n
_r = n / (2.0 * n + 2.0)

def H0(t):
    return _H0 * (1.0 - 0.5 * np.sin(np.pi * t / _T))

def L(t):
    return _L * (1.0 - 0.75 * np.sin(np.pi * t / _T))

def _psi(t,x):   # vectorized over x
    psi = np.zeros(np.shape(x))
    y = np.abs(x[np.abs(x)<L(t)]) / L(t)
    psi[np.abs(x)<L(t)] = (n+1.0) * y - 1.0 + n * (1.0 - y)**_q - n * y**_q
    return psi

def s(t,x):   # vectorized over x
    _C = H0(t) / (n - 1.0)**_r
    return _C * _psi(t,x)**_r

def sx(t,x):  # vectorized over x
    _C = _r * H0(t) / (n - 1.0)**_r
    gamma = np.zeros(np.shape(x))
    gamma[np.abs(x)<L(t)] = _psi(t,x[np.abs(x)<L(t)])**(_r-1.0)
    omega = np.zeros(np.shape(x))
    y = np.abs(x[np.abs(x)<L(t)]) / L(t)
    z = np.sign(x[np.abs(x)<L(t)]) / L(t)
    omega[np.abs(x)<L(t)] = (n+1.0) * z - (n+1.0) * (1.0 - y)**(1.0/n) * z \
                            - (n+1.0) * y**(1.0/n)  * z
    return _C * gamma * omega

def u(t,x):  # vectorized over x
    _gam = 2.0 * A * (rho * g)**n / (n+1)
    return - _gam * s(t,x)**(n+1.0) * sx(t,x)**n   # assumes n odd integer

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    print('plotting glacier at t=0, %.0f, %.0f years ...' \
          % (0.25*_T/secpera,0.5*_T/secpera))
    xx = np.linspace(-1.1*_L,1.1*_L,801)
    plt.subplot(3,1,1)
    for tt in [0.0, 0.25*_T, 0.5*_T]:
        plt.plot(xx/1.0e3,s(tt,xx),
                 label='%.0f years' % (tt/secpera))
    plt.ylabel('elevation s  (m)')
    plt.legend()
    plt.subplot(3,1,2)
    for tt in [0.0, 0.25*_T, 0.5*_T]:
        plt.plot(xx/1.0e3,sx(tt,xx),'.',ms=2.0)
    plt.ylabel('surface slope s_x')
    plt.subplot(3,1,3)
    for tt in [0.0, 0.25*_T, 0.5*_T]:
        plt.plot(xx/1.0e3,u(tt,xx)*secpera)
    plt.ylabel('surface velocity u  (m/a)')
    plt.xlabel('x  (km)')
    plt.show()
    # TODO: compute ds/dt, a
