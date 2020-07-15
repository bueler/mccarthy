from firedrake import *

mesh = UnitSquareMesh(4,4)

X = FunctionSpace(mesh, 'CG', 1)
a = Function(X)
File('foo.pvd').write(a)

V = VectorFunctionSpace(mesh, 'CG', 2)
W = FunctionSpace(mesh, 'CG', 1)
Z = V * W
up = Function(Z)
#u,p = split(up)
u,p = up.split()
File('bar.pvd').write(u,p)

