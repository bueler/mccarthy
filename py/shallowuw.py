'''
Compute and plot the surface values of the horizontal
and vertical velocity, (u,w), using the non-sliding
shallow ice approximation (SIA).
'''

import numpy as np
import matplotlib.pyplot as plt
secpera = 31556925.0

# major parameters
g = 9.81               # acceleration of gravity (m s-2)
rhoi = 910.0           # density of ice (kg m-3)
n = 3                  # Glen exponent
A = 3.1689e-24         # EISMINT I value of ice softness (Pa-3 s-1)
L = 20.0e3             # length of domain (m)
H0 = 500.0             # reference ice thickness (m)
alpha = 0.1            # slope of bed (radians)
J = 40                 # number of subintervals (in x direction)
outdir = 'output/'     # directory name for image output

def s(x, J=J, amplitude=100.0, wavelength=4.0e3):
    '''compute s(x) with sine wave'''
    return H0 + amplitude * np.sin(2.0 * np.pi * x / wavelength)

def compute_velocity(x, s, b):
    '''compute surface values of velocity u(x), w(x) from
    s and b using SIA formula'''
    dx = x[1] - x[0]
    dsdx = (s1d[1:] - s1d[:-1]) / dx
    u = - 0.5 * A * (rhoi * g * dsdx)**n * (s - b)**(n+1)
    FIXME
    return u, w

def mkoutdir(dirname):
    import os
    try:
        os.mkdir(dirname)
    except FileExistsError:
        pass

def plot_velocity():
    FIXME
    # save plot to .png
    xplot = np.linspace(0.0, L, 1001)
    plt.figure(2, figsize=(12,5))
    plt.plot(xplot, b(xplot))
    plt.axis([0, L, 0, 1.5 * b(0.0)])
    plt.xlabel('x (m)')
    plt.ylabel('elevation (m)')
    plt.plot(x, s, "-o")
    plt.text(1000.0, 700.0, 'steady state', name='DejaVu Sans Mono')
    mkoutdir(outdir)
    print(f'  writing to image file {outdir}steady.png')
    plt.savefig(f'{outdir}steady.png')

FIXME
u, w = compute_velocity()
plot_velocit(s, u, w)
