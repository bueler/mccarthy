import numpy as np
from scipy import linalg

m, n = 3, 2
A = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
b = np.array([3.141, 3.0, 2.718])
print(A)
print(b)

x, res, rnk, s = linalg.lstsq(A, b)  # res, rnk, s give extra info
print(x)

r = np.dot(A,x) - b
print(r)
print(linalg.norm(r)/linalg.norm(b))
