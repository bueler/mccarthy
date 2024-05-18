# generate coefficients in b(x) polynomial:
#   b(x) = a2 * (x-z)**2 + a3 * (x-z)**3 + a4 * (x-z)**4 + a5 * (x-z)**5
# where z = L/2, so that b(z)=b'(z)=0 and
#     b0 = b(0)
#    -m0 = b'(0)
#   b0/3 = b(L)
#      0 = b'(L)
# in fact a5=0, so b(x) is quartic

import numpy as np
import matplotlib.pyplot as plt

L, b0, m0 = 20.0e3, 500.0, 0.1

z = L/2
Amat = np.array([[z**2,  -z**3,    z**4,  -z**5],
                 [-2*z, 3*z**2, -4*z**3, 5*z**4],
                 [z**2,   z**3,    z**4,   z**5],
                 [ 2*z, 3*z**2,  4*z**3, 5*z**4]])
brhs = np.array([b0, -m0, b0/3, 0.0])
a2, a3, a4, a5 = tuple(np.linalg.solve(Amat,brhs))
assert abs(a5 * z**3 / a2) < 1.0e-14
a5 = 0.0

# result: this block computes b(x)
print(f'a2={a2:.3e}, a3={a3:.3e}, a4={a4:.3e}')
def b(x):
    x = np.array(x)
    z = L / 2
    return (x-z)**2 * (a2 + a3 * (x-z) + a4 * (x-z)**2)

# visualize to confirm
x = np.linspace(-1e3, L+1e3, 1001)
plt.plot(x, b(x))
plt.plot([0,L], b([0,L]), 'b+')
plt.show()

