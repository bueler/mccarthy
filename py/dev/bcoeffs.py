# generate coefficients in b(x) polynomial:
#   b(x) = c2 * (x-z)**2 + c3 * (x-z)**3 + c4 * (x-z)**4 + c5 * (x-z)**5
# where z = L/2, so that b(z)=b'(z)=0 and
#     b0 = b(0)
#    -m0 = b'(0)
#   b0/3 = b(L)
#      0 = b'(L)
# in fact c5=0, so b(x) is quartic

import numpy as np
import matplotlib.pyplot as plt

L, b0, m0 = 20.0e3, 500.0, 0.1

# solve linear system and print coefficients
z = L/2
Amat = np.array([[z**2,  -z**3,    z**4,  -z**5],
                 [-2*z, 3*z**2, -4*z**3, 5*z**4],
                 [z**2,   z**3,    z**4,   z**5],
                 [ 2*z, 3*z**2,  4*z**3, 5*z**4]])
brhs = np.array([b0, -m0, b0/3, 0.0])
c2, c3, c4, c5 = tuple(np.linalg.solve(Amat, brhs))
assert abs(c5 * z**3 / c2) < 1.0e-14
c5 = 0.0
print(f'c2, c3, c4 = {c2:.3e}, {c3:.3e}, {c4:.3e}')

# result
def b(x):
    x = np.array(x)
    z, c2, c3, c4 = L / 2, 4.167e-06, -1.667e-10, -8.333e-15
    return (x-z)**2 * (c2 + c3 * (x-z) + c4 * (x-z)**2)

# visualize to confirm
x = np.linspace(-1e3, L+1e3, 1001)
plt.plot(x, b(x))
plt.plot([0,L], b([0,L]), 'b+')
plt.show()
