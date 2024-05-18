import numpy as np
import matplotlib.pyplot as plt

L = 20.0e3
J = 20
secpera = 31556925.0
dt = 3.0 * secpera


def a(x):
    return 1.0e-7 - 2.0e-7 * x / L


def b(x):
    z, c2, c3, c4 = L / 2, 4.167e-06, -1.667e-10, -8.333e-15
    return (x - z) ** 2 * (c2 + c3 * (x - z) + c4 * (x - z) ** 2)


def u(x):
    return 2.0e-7 + 0.0 * x


x = np.linspace(0.0, L, J + 1)
dx = x[1] - x[0]
xplot = np.linspace(0.0, L, 1001)

# evolution loop
s = b(0.0) - 200.0 * x / L + 50.0 * np.sin(10.0 * x / L)
plt.plot(xplot, b(xplot))
plt.plot(x, s, "-o")
for k in range(4):
    snew = s.copy()
    snew[1:] += dt * (a(x[:-1]) - u(x[:-1]) * (s[1:] - s[:-1]) / dx)
    s = snew
    plt.plot(x, s, "-o")

plt.show()
