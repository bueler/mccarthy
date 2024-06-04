'''
Compute and plot the surface values of the horizontal and vertical
velocity using the non-sliding shallow ice approximation (SIA).
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
alpha = 0.1            # slope of bed
wave_amp = 100.0       # amplitude of surface waves
wave_len = 4.0e3       # wavelength of surface waves
J = 40                 # number of subintervals (in x direction)
outdir = 'output/'     # directory name for image output

# derived constant from SIA formulas
C = (2.0/(n+1)) * A * (rhoi * g)**n

def b(x):
    'Compute b(x), a slope.'
    return - alpha * x

def s(x):
    'Compute s(x), a sine wave.'
    return b(x) + H0 + wave_amp * np.sin(2.0 * np.pi * x / wave_len)

def u(x, b, s):
    '''Compute surface value of horizontal velocity u = u_j at
    interior regular grid locations.'''
    dx = x[1] - x[0]
    dsdx = (s[1:] - s[:-1]) / dx
    ult = - C * abs(dsdx[:-1])**(n-1) * dsdx[:-1] * (s[:-1] - b[:-1])**(n+1)
    urt = - C * abs(dsdx[1:] )**(n-1) * dsdx[1:]  * (s[1:]  - b[1:] )**(n+1)
    return (ult + urt) / 2.0

def vert_int_u_stag(j, x, b, s, zb, zs):
    '''Vertically integrate the horizontal velocity at a
    staggered location jj = j+1/2 from level zb up to level zs.
    Inputs are a horizontal index j, three vectors x, b, s
    (functions of x), and two scalars zb, zs.

    This most-technical code integrates ujj(z), which is
    defined piecewise as follows.  For a given j the staggered-
    location values of b, s, and surface slope are computed
    (namely: bjj, sjj, dsjj).  SIA formulas convert these to
    constant values alf, bet in the formula for the horizontal
    velocity:
                / 0
      ujj(z) = |  alf * (bet - (sjj - z)^(n+1))
                \ alf * bet
    for vertical ranges z <= bjj, bjj < z < sjj, sjj <= z
    respectively.  Note that ujj = 0 below the bed and that
    it is continued above the surface by its surface value.

    The returned value is the exact integral
      I = int_zb^zs ujj(z) dz.  Six cases are needed to
    compute this integral.
    '''
    assert zb <= zs
    dx = x[1] - x[0]
    bjj = (b[j] + b[j+1]) / 2.0
    sjj = (s[j] + s[j+1]) / 2.0
    dsjj = (s[j+1] - s[j]) / dx
    alf = - 0.5 * A * (rhoi * g)**n * abs(dsjj)**(n-1) * dsjj
    bet = (sjj - bjj)**(n+1)
    if zb <= bjj
        if zs <= bjj:
            I = 0.0
        else if zs < sjj:
            I = FIXME
        else:
            I = FIXME + alf * bet * (zs - sjj)
    else if zb < sjj:
        if zs < sjj:
            I = FIXME
        else:
            I = FIXME + alf * bet * (zs - sjj)
    else:
        I = alf * bet * (zs - zb)
    return I

def w(x, b, s):
    '''Compute surface value of vertical velocity w = w_j at
    interior regular grid locations.'''
    dx = x[1] - x[0]
    Ilt = np.zeros(shape=(J-1,))
    Irt = np.zeros(shape=(J-1,))
    for j in range(J-1):
        Ilt[j] = vert_int_u_stag(j,   x, b, s, b[j+1], s[j+1])
        Irt[j] = vert_int_u_stag(j+1, x, b, s, b[j+1], s[j+1])
    return - (Ilt - Irt) / dx

def mkoutdir(dirname):
    import os
    try:
        os.mkdir(dirname)
    except FileExistsError:
        pass

def plot_velocity():
    x = np.linspace(0.0, L, J+1)
    xint = x[1:-1]
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(xint, u(x, b(x), s(x)))
    ax2.plot(xint, w(x, b(x), s(x)))
    ax2.xlabel('x (m)')
    #plt.axis([0, L, 0, 1.5 * b(0.0)])
    #plt.ylabel('elevation (m)')
    #plt.text(1000.0, 700.0, 'steady state', name='DejaVu Sans Mono')
    mkoutdir(outdir)
    print(f'  writing to image file {outdir}uw.png')
    plt.savefig(f'{outdir}uw.png')

plot_velocity()
