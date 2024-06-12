from firedrake import *

# define mesh and function spaces on [0,L]^2
L, N = 50.0e3, 100
mesh = RectangleMesh(N, N, L, L)
P1 = FunctionSpace(mesh, "CG", 1)
P1vector = VectorFunctionSpace(mesh, "CG", 1, dim=3)

# define data a, b, u
secpera = 31556926.0

def surfacemassbalance(xy):
    lapserate = 0.2 / (1.0e3 * secpera)
    maxsmb = 0.2 / secpera
    center = as_vector([15e3, 40e3])
    R = sqrt(dot(xy - center, xy - center))
    return maxsmb - lapserate * R

def bedelevation(xy):
    return 600.0 - 300.0 * (L - xy[1] + xy[0])/L \
           + 200.0 * ((L - xy[1] - xy[0])/L)**2

def surfacevelocity(xy):
    U0 = 300.0 / secpera
    return as_vector([Constant(U0), -Constant(U0), 0.0])

# define (regularized) weak form problem
xy = SpatialCoordinate(mesh)
a = Function(P1, name='a(x,y)').interpolate(surfacemassbalance(xy))
b = Function(P1, name='b(x,y)').interpolate(bedelevation(xy))
u = Function(P1vector, name='u(x,y)').interpolate(surfacevelocity(xy))
s = Function(P1, name='s(x,y)').interpolate(b)  # initialize to no ice
ns = as_vector([-s.dx(0), -s.dx(1), 1.0])
v = TestFunction(P1)
epsreg = Constant(0.001)
F = epsreg * inner(grad(s), grad(v)) * dx \
    + (-inner(u, ns) - a) * v * dx
bcs = DirichletBC(P1, b, "on_boundary")
problem = NonlinearVariationalProblem(F, s, bcs)

# solve variational inequality (VI) problem
sp = {"snes_type": "vinewtonrsls",
      "snes_converged_reason": None,
      "snes_linesearch_type": "basic",
      "ksp_type": "preonly",
      "pc_type": "lu",
      "pc_factor_mat_solver_type": "mumps"}
nvs = NonlinearVariationalSolver(problem,
                                    solver_parameters=sp)
upper = Function(P1).interpolate(Constant(10000.0))
nvs.solve(bounds=(b, upper))

# Paraview-readable output
from firedrake.output import VTKFile
H = Function(P1, name='H(x,y)').interpolate(s - b)  # ice thickness
VTKFile('output/surface2d.pvd').write(a,b,u,s,H)

# image output (greyscale contour map)
import numpy as np
import matplotlib.pyplot as plt
from firedrake.pyplot import tricontour
fig, axes = plt.subplots()
tricontourf(H, axes=axes, levels=[0.0, 1.0], cmap="Greys")
contours = tricontour(s, axes=axes, levels=15, cmap="magma")
fig.colorbar(contours, label='s (m)')
axes.set_aspect("equal")
tickskm = np.linspace(0,L/1000.0,6)
tickstrs = ['0','10','20','30','40','50']  # for L = 50 km
axes.set_xticks(tickskm * 1000.0, tickstrs)
plt.xlabel('x (km)')
axes.set_yticks(tickskm * 1000.0, tickstrs)
plt.ylabel('y (km)')
plt.savefig('output/surface2d.png', bbox_inches='tight')