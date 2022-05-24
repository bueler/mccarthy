import numpy as np
from scipy import linalg

m, n = 3, 2
A = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
print(A)

U, s, Vt = linalg.svd(A, full_matrices=False)
Sigma = np.zeros((n,n))
np.fill_diagonal(Sigma, s)
print(Sigma)
print(linalg.norm(A - np.dot(U, np.dot(Sigma,Vt))))
