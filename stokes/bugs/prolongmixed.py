from firedrake import *

mesh = UnitSquareMesh(4,4)
hierarchy = MeshHierarchy(mesh, 1)

V0 = VectorFunctionSpace(hierarchy[0], 'CG', 2)
W0 = FunctionSpace(hierarchy[0], 'CG', 1)
Z0 = V0 * W0
up0 = Function(Z0)
#u0,p0 = split(up0)
u0,p0 = up0.split()

V1 = VectorFunctionSpace(hierarchy[1], 'CG', 2)
W1 = FunctionSpace(hierarchy[1], 'CG', 1)
Z1 = V1 * W1
up1 = Function(Z1)
#u1,p1 = split(up1)
u1,p1 = up1.split()

prolong(u0,u1)
prolong(p0,p1)

