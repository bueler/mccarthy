from firedrake import *

# FIXME make L smaller for a valley glacier, e.g. 50km

# define mesh and function spaces on [0,L]^2
L, N = 1000.0e3, 50
mesh = RectangleMesh(N, N, L, L)
V = FunctionSpace(mesh, "CG", 1)
W = VectorFunctionSpace(mesh, "CG", 1, dim=3)

# define data a, b, u
secpera = 31556926.0

def surfacemassbalance(xy):
    lapserate = 0.0004 / (1.0e3 * secpera)
    maxsmb = 0.2 / secpera
    center = as_vector([-100e3, 900e3])
    R = sqrt(dot(xy - center, xy - center))
    return maxsmb - lapserate * R

def bedelevation(xy):
    return 600.0 - 300.0 * (L - xy[1] + xy[0])/L \
           + 200.0 * ((L - xy[1] - xy[0])/L)**2

def surfacevelocity(xy):
    U0 = 100.0 / secpera
    return as_vector([Constant(U0), -Constant(U0), 0.0])

# define (regularized) weak form problem
xy = SpatialCoordinate(mesh)
a = Function(V, name='a(x,y)').interpolate(surfacemassbalance(xy))
b = Function(V, name='b(x,y)').interpolate(bedelevation(xy))
u = Function(W, name='u(x,y)').interpolate(surfacevelocity(xy))
s = Function(V, name='s(x,y)').interpolate(b)  # initialize to no ice
ns = as_vector([-s.dx(0), -s.dx(1), 1.0])
v = TestFunction(V)
epsreg = Constant(0.01)
F = epsreg * inner(grad(s), grad(v)) * dx \
    + (-inner(u, ns) - a) * v * dx
bcs = DirichletBC(V, b, "on_boundary")
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
upper = Function(V).interpolate(Constant(10000.0))
nvs.solve(bounds=(b, upper))

from firedrake.output import VTKFile
VTKFile('mountainglacier.pvd').write(a,b,u,s)
