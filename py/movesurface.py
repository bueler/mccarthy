import numpy as np

L = 20.0e3
J = 20
secpera = 31556925.0
dt = 1.0

x = np.linspace(0.0, L, J+1)
xplot = np.linspace(0.0, L, 1001)

b = 0.0 * x # FIXME design it
s = 0.0 * x # FIXME data for initial surface

#snew = s + 