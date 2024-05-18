# generate coefficients in b(x) polynomial:
#   b(x) = a2 * (x-z)**2 + a3 * (x-z)**3 + a4 * (x-z)**4 + a5 * (x-z)**5
# where z = L/2

# FIXME progress, but redesign

import numpy as np
import matplotlib.pyplot as plt

L, b0, m0 = 20.0e3, 500.0, 0.1

z, w = L/2, L/4
Amat = np.array([[z**2,  -z**3,    z**4,  -z**5],
                 [-2*z, 3*z**2, -4*z**3, 5*z**4],
                 [w**2,   w**3,    w**4,   w**5],
                 [2*w,  3*w**2,  4*w**3, 5*w**4]])
brhs = np.array([b0, -m0, b0/2, 0.0])
a2, a3, a4, a5 = tuple(np.linalg.solve(Amat,brhs))

x = np.linspace(0.0, L, 1001)
b = a2  * (x-z)**2 + a3 * (x-z)**3 + a4 * (x-z)**4 + a5 * (x-z)**5
plt.plot(x, b)
plt.show()

