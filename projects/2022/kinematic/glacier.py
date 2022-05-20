#!/usr/bin/env python3
# (C) 2022 Ed Bueler

import numpy as np

# t in seconds, x in meters, s(t,x) = H(t,x) in meters, u(t,x) in m/s
secpera = 31556926.0
n = 3.0
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
    y = np.minimum(1.0, np.abs(x) / L(t))
    return (n+1.0) * y - 1.0 + n * (1.0 - y)**_q - n * y**_q

def s(t,x):   # vectorized over x
    _C = H0(t) / (n - 1.0)**_r
    return _C * _psi(t,x)**_r

def _sx(t,x):  # vectorized over x
    _C = _r * H0(t) / (n - 1.0)**_r
    y = np.abs(x) / L(t)
    #FIXME return _C * _psi(t,x) *
    return 0.0

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    print('plotting glacier at t=0, %.0f, %.0f years ...' \
          % (0.25*_T/secpera,0.5*_T/secpera))
    xx = np.linspace(-1.1*_L,1.1*_L,401)
    for tt in [0.0, 0.25*_T, 0.5*_T]:
        ss = s(tt,xx)
        plt.plot(xx/1.0e3,ss,label='%.0f years' % (tt/secpera))
    plt.xlabel('x  (km)')
    plt.ylabel('s  (m)')
    plt.legend()
    plt.show()
    # TODO: compute s_x, u, ds/dt, a
