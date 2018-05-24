# module for mesh-related computations/actions

from firedrake import Function, FunctionSpace, VectorFunctionSpace, \
                      DirichletBC, Constant, SpatialCoordinate, \
                      TrialFunction, TestFunction, grad, inner, dx, solve
from scipy.interpolate import interp1d

# return integer-valued fields on the mesh which give the process rank for
# element ownership (piecewise constant; discontinuous) and vertex/node
# ownership (integer-valued on vertices; piecewise linear on mesh)
def getranks(mesh):
    element_rank = Function(FunctionSpace(mesh,'DG',0))
    element_rank.dat.data[:] = mesh.comm.rank
    element_rank.rename('element_rank')
    vertex_rank = Function(FunctionSpace(mesh,'CG',1))
    vertex_rank.dat.data[:] = mesh.comm.rank
    vertex_rank.rename('vertex_rank')
    return (element_rank,vertex_rank)

# return linear-interpolated surface elevation function  z = h(x)  for current mesh
# notes: 1)  bc.nodes  gives indices to mesh; is a 1D dtype=int32 numpy array
#        2)  bc.nodes  includes the halo nodes so need to use .dat.data_with_halos
#        3)  f = Function(P1);  bc.apply(f)  would give indicator function
#        4)  default for interp1d() is 'linear', which is what we want
#        5)  by turning off bounds_error, and allowing extrapolation,
#            there is parallel functionality for this function
def getsurfaceelevationfunction_halos(mesh,top_id):
    P1 = FunctionSpace(mesh, "CG", 1)
    bc = DirichletBC(P1, 1.0, top_id)
    x,z = SpatialCoordinate(mesh)
    xh_halos = Function(P1).interpolate(x).dat.data_with_halos[bc.nodes]  # 1D numpy array
    zh_halos = Function(P1).interpolate(z).dat.data_with_halos[bc.nodes]
    return interp1d(xh_halos,zh_halos,copy=False,bounds_error=False,fill_value='extrapolate')

# compute vertical mesh displacement, given the surface value of it, by setting
# up and solving a Dirichlet (Laplace equation) problem in parallel
def solvevdisp(mesh,bdryids,deltah):
    P1 = FunctionSpace(mesh, "CG", 1)
    r = TrialFunction(P1)
    s = TestFunction(P1)
    a = inner(grad(r), grad(s)) * dx   # note natural b.c. on outflow
    L = inner(Constant(0.0), s) * dx
    # WARNING: top must go *first* so closed top gets zero; is this documented behavior?
    bcs = [ DirichletBC(P1, deltah, bdryids['top']),
            DirichletBC(P1, Constant(0.0), (bdryids['base'],bdryids['inflow'])) ]
    rsoln = Function(P1)
    solve(a == L, rsoln, bcs=bcs, options_prefix='t', solver_parameters={})
    return rsoln

# following functions are serial-only, and for generating plots

# helper function
def getxsurface(mesh,top_id):
    P1 = FunctionSpace(mesh, "CG", 1)
    bc = DirichletBC(P1, 1.0, top_id)
    x,_ = SpatialCoordinate(mesh)
    xh = Function(P1).interpolate(x).dat.data_ro[bc.nodes]
    return (xh, bc)

# return linear-interpolated surface vertical displacement function  r = (delta h)(x)
def getsurfacevdispfunction(mesh,top_id,r):
    xh,bc = getxsurface(mesh,top_id)
    return interp1d(xh,r.dat.data_ro[bc.nodes],copy=False)

# return linear-interpolated surface velocity functions  u(x), w(x)
def getsurfacevelocityfunction(mesh,top_id,Z,u):
    xh,bc = getxsurface(mesh,top_id)
    P1V = VectorFunctionSpace(mesh, "CG", 1)
    uP1 = Function(P1V).interpolate(u)
    ufcn = interp1d(xh,uP1.dat.data_ro[bc.nodes,0],copy=False)
    wfcn = interp1d(xh,uP1.dat.data_ro[bc.nodes,1],copy=False)
    return (ufcn,wfcn)

