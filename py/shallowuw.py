'''
Compute and plot the surface values of the horizontal and vertical
velocity using the non-sliding shallow ice approximation (SIA).
Techniques emphasize robustness, e.g. when ice thickness goes to
zero.
'''

import numpy as np
import matplotlib.pyplot as plt
secpera = 31556925.0

wave = True            # set to False to compare slab-on-slope formula

# instance parameters
L = 20.0e3             # length of domain (m)
H0 = 400.0             # reference ice thickness (m)
alpha = -0.02          # slope of bed
if wave:
    wave_amp = 20.0    # amplitude of surface waves
else:
    wave_amp = 0.0
wave_len = 8.0e3       # wavelength of surface waves
J = 80                 # number of subintervals (in x direction)
outdir = 'output/'     # directory name for image output

# physical parameters
g = 9.81               # acceleration of gravity (m s-2)
rhoi = 910.0           # density of ice (kg m-3)
n = 3                  # Glen exponent
A = 3.1689e-24         # EISMINT I value of ice softness (Pa-3 s-1)

# derived constant from SIA formulas
C = (2.0 / (n+1)) * A * (rhoi * g)**n

def b(x):
    'Compute b(x), a slope.'
    return alpha * (x - L)

def s(x):
    'Compute s(x), a sine wave, but satisfying s(x) >= b(x).'
    stmp = H0 + wave_amp * np.sin(2.0 * np.pi * x / wave_len)
    return b(x) + np.maximum(stmp, 0.0)

def u(x, b, s):
    '''Compute surface value of horizontal velocity u = u_j at
    interior regular grid locations.'''
    assert np.all(b <= s)
    dx = x[1] - x[0]
    dsdx = (s[1:] - s[:-1]) / dx
    Hstag = ((s[:-1] - b[:-1]) + (s[1:] - b[1:])) / 2.0
    ult = - C * abs(dsdx[:-1])**(n-1) * dsdx[:-1] * Hstag[:-1]**(n+1)
    urt = - C * abs(dsdx[1:])**(n-1)  * dsdx[1:]  * Hstag[1:]**(n+1)
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
    assert bjj <= sjj
    dsjj = (s[j+1] - s[j]) / dx
    alf = - 0.5 * A * (rhoi * g)**n * abs(dsjj)**(n-1) * dsjj
    bet = (sjj - bjj)**(n+1)
    if zb <= bjj:
        if zs <= bjj:
            I = 0.0
        elif zs < sjj:
            I = alf * bet * (zs - bjj) \
                - (alf / (n+2)) * ((sjj - zs)**(n+2) - (sjj - bjj)**(n+2))
        else:
            I = alf * bet * (sjj - bjj) \
                - (alf / (n+2)) * (0.0 - (sjj - bjj)**(n+2)) \
                + alf * bet * (zs - sjj)
    elif zb < sjj:
        if zs < sjj:
            I = alf * bet * (zs - zb) \
                - (alf / (n+2)) * ((sjj - zs)**(n+2) - (sjj - zb)**(n+2))
        else:
            I = alf * bet * (sjj - zb) \
                - (alf / (n+2)) * (0.0 - (sjj - zb)**(n+2)) \
                + alf * bet * (zs - sjj)
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

def compute_velocity():
    x = np.linspace(0.0, L, J+1)
    uu = u(x, b(x), s(x)) * secpera
    ww = w(x, b(x), s(x)) * secpera
    print(f'  range for u:  min = {uu.min():9.4f},  max = {uu.max():9.4f}  (m/a)')
    print(f'  range for w:  min = {ww.min():9.4f},  max = {ww.max():9.4f}  (m/a)')
    if not wave:
        # in slab-on-slope case, compare horizontal surface velocity
        u_slab = - (2.0 / (n+1)) * A * (rhoi * g)**n * abs(alpha)**2 * alpha * H0**(n+1)
        print(f'  COMPARE: slab-on-slope says u={u_slab * secpera:.4f} (m/a)')
    return x, uu, ww

def plot_velocity(x, u, w):
    xint = x[1:-1]
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(x / 1.0e3, s(x), '.-', color='C1', ms=4.0, label='s')
    ax1.plot(x / 1.0e3, b(x), '-', color='C0', label='b')
    ax1.legend(loc='lower left')
    ax1.set_xticklabels([])
    ax2.set_ylabel('elevation (m)')
    ax2.plot(xint / 1.0e3, u, 'o', ms=6.0, mfc='w', label='u')
    ax2.plot(xint / 1.0e3, w, '+', ms=8.0, label='w')
    ax2.legend(loc='upper right')
    ax2.set_xticks(np.linspace(0.0,20.0,5))
    ax2.set_xticklabels(['0','5','10','15','20'])
    ax2.set_ylabel('velocity (m/a)')
    ax2.grid(visible=True)
    plt.xlabel('x (km)')
    mkoutdir(outdir)
    print(f'  writing to image file {outdir}uw.png')
    plt.savefig(f'{outdir}uw.png')

x, u, w = compute_velocity()
plot_velocity(x, u, w)
