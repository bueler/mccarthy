# module for motion of mesh from surface kinematical equation

import firedrake as fd
from surfaceutils import getsurfaceelevation

# use surface kinematical equation to compute mesh vertical displacement field
#   computes vertical mesh displacement given its surface value
#   sets-up and solves a Dirichlet (Laplace equation) problem in parallel
#   note solver prefix is -vd_
#   POSSIBILITY TO IMPROVE?:  make the Laplacian isotropic
#   POSSIBILITY: use better solver on Laplace equation
#   POSSIBILITY: add in climatic mass balance a(x) here
#                h_t = a - u[0] h_x + u[1], but currently uses a = Constant(0.0)
def surfacekinematical(mesh,bdryids,u,dt):
    # get deltah = change in surface elevation
    h = getsurfaceelevation(mesh,bdryids['top'])
    x,z = fd.SpatialCoordinate(mesh)
    P1 = fd.FunctionSpace(mesh,'CG',1)
    xval = fd.Function(P1).interpolate(x)
    zval = fd.Function(P1).interpolate(z)
    phi = fd.Function(P1)
    phi.dat.data[:] = zval.dat.data_ro - h(xval.dat.data_ro)
    deltah = fd.Function(P1).interpolate( dt * (fd.Constant(0.0) + fd.inner(fd.grad(phi),u)) )
    # solve mesh displacement problem by Laplace equation
    r = fd.TrialFunction(P1)
    s = fd.TestFunction(P1)
    a = fd.inner(fd.grad(r), fd.grad(s)) * fd.dx   # note natural b.c. on outflow
    L = fd.inner(fd.Constant(0.0), s) * fd.dx
    # WARNING: top must go *first* so closed top gets zero; is this documented behavior?
    bcs = [ fd.DirichletBC(P1, deltah, bdryids['top']),
            fd.DirichletBC(P1, fd.Constant(0.0), (bdryids['base'],bdryids['inflow'])) ]
    rsoln = fd.Function(P1)
    fd.solve(a == L, rsoln, bcs=bcs, options_prefix='vd', solver_parameters={})
    rsoln.rename('vertical_displacement')
    return rsoln

# move mesh by applying vertical displacement field r
def movemesh(mesh,r,bmin):
    x,z = fd.SpatialCoordinate(mesh)
    Vc = mesh.coordinates.function_space()
    f = fd.Function(Vc).interpolate(fd.as_vector([x, z + r]))
    mesh.coordinates.assign(f)
    unstable = any(f.dat.data_ro[:,1] < bmin - 1.0)
    return mesh,unstable

