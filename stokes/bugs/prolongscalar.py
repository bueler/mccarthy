from firedrake import *

mesh = UnitSquareMesh(4,4)
hierarchy = MeshHierarchy(mesh, 1)

W0 = FunctionSpace(hierarchy[0], 'CG', 1)
p0 = Function(W0)

W1 = FunctionSpace(hierarchy[1], 'CG', 1)
p1 = Function(W1)

prolong(p0,p1)

