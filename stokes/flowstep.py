#!/usr/bin/env python3
# (C) 2018 Ed Bueler

# Solve glacier bedrock-step Glen-Stokes problem.  See
# mccarthy/stokes/README.md for usage.

import argparse
secpera = 31556926.0       # seconds per year

# process options
mixchoices = ['P2P1','P3P2','P2P0','CRP0','P1P0']
parser = argparse.ArgumentParser(\
    description='Solve glacier Glen-Stokes problem on bedrock-step geometry.')
parser.add_argument('-alpha', type=float, default=0.1, metavar='X',
                    help='downward slope of bed as angle in radians (default = 0.1)')
parser.add_argument('-Dtyp', type=float, default=2.0, metavar='X',
                    help='regularize viscosity using "+(eps Dtyp)^2" (default = 2 a-1)')
                    # e.g. (800 m a-1) / 400 m = 2 a-1
parser.add_argument('-elements', metavar='X', default='P2P1',
                    choices=mixchoices,
                    help='stable mixed finite elements from: %s (default=P2P1)' \
                         % (','.join(mixchoices)) )
parser.add_argument('-eps', type=float, default=0.01, metavar='X',
                    help='regularize viscosity using "+(eps Dtyp)^2" (default = 0.01)')
parser.add_argument('inname', metavar='INNAME',
                    help='input file name ending with .msh')
parser.add_argument('-n_glen', type=float, default=3.0, metavar='X',
                    help='Glen flow law exponent (default = 3.0)')
parser.add_argument('-o', metavar='OUTNAME', type=str, default='',
                    help='output file name ending with .pvd (default = INNAME-.msh+.pvd)')
args, unknown = parser.parse_known_args()
if len(args.o) > 0:
    outname = args.o
else:
    outname = '.'.join(args.inname.split('.')[:-1]) + '.pvd'  # strip .msh and replace with .pvd
n_glen = args.n_glen  # used pretty often
Dtyp = args.Dtyp / secpera

from firedrake import *

# glacier physical constants
g = 9.81                   # m s-2
rho = 910.0                # kg m-3
A3 = 3.1689e-24            # Pa-3 s-1; EISMINT I value of ice softness
B3 = A3**(-1.0/3.0)        # Pa s(1/3);  ice hardness

# input mesh and define geometry
print('reading mesh from %s ...' % args.inname)
mesh = Mesh(args.inname)
print('mesh has %d vertices and %d elements' \
      % (mesh.num_vertices(),mesh.num_cells()))

# numbering of parts of boundary *must match generation script genstepmesh.py*
outflow_id = 41
top_id = 42  # unused below
inflow_id = 43
base_id = 44

# extract mesh geometry making these assumptions (tolerance=1cm):
#   * max x-coordinate is total length L
#   * min z-coordinate at x=0 is bs
#   * max z-coordinate at x=0 is hsurfin
#   * max z-coordinate at x=L is Hout
xarray = mesh.coordinates.dat.data[:,0]
zarray = mesh.coordinates.dat.data[:,1]
L = max(xarray)
bs = min(zarray[xarray<0.01])
hsurfin = max(zarray[xarray<0.01])
Hin = hsurfin - bs
Hout = max(zarray[xarray>L-0.01])
print('mesh geometry [m]: L = %.3f, bs = %.3f, Hin = %.3f, Hout = %.3f' \
      %(L,bs,Hin,Hout))
isslab = False
if abs(bs) < 0.01 and abs(Hin - Hout) < 0.01:
    print('  ... detected slab geometry case ...')
    isslab = True
print('using bed slope angle alpha = %.6f' % args.alpha)

# determine B_n so that slab-on-slope solutions give surface velocity that is n-independent
Bn = (4.0/(n_glen+1.0))**(1.0/n_glen) \
     * (rho*g*sin(args.alpha)*Hin)**((n_glen-3.0)/n_glen) * B3**(3.0/n_glen)

# define mixed finite element space
mixFE = {'P2P1'  : (VectorFunctionSpace(mesh, "CG", 2), # Taylor-Hood
                    FunctionSpace(mesh, "CG", 1)),
         'P3P2'  : (VectorFunctionSpace(mesh, "CG", 3),
                    FunctionSpace(mesh, "CG", 2)),
         'P2P0'  : (VectorFunctionSpace(mesh, "CG", 2),
                    FunctionSpace(mesh, "DG", 0)),
         'CRP0'  : (VectorFunctionSpace(mesh, "CR", 1),
                    FunctionSpace(mesh, "DG", 0)),
         'P1P0'  : (VectorFunctionSpace(mesh, "CG", 1), # interesting but not
                    FunctionSpace(mesh, "DG", 0))       #   recommended
        }
(V,W) = mixFE[args.elements]
Z = V * W
v,q = TestFunctions(Z)

# define body force
f_body = Constant((g * rho * sin(args.alpha), - g * rho * cos(args.alpha)))

# Dirichlet boundary condition applied on base
noslip = Constant((0.0, 0.0))

# right side outflow: apply hydrostatic normal force; nonhomogeneous Neumann
x,z = SpatialCoordinate(mesh)
outflow_sigma = as_vector([- rho * g * cos(args.alpha) * (Hout - z),
                           rho * g * sin(args.alpha) * (Hout - z)])

# put solution here
up = Function(Z)       # *not* TrialFunctions(Z)
u,p = split(up)        # up.split() not equivalent here?

def D(w):
    return 0.5 * (grad(w) + grad(w).T)

# define the nonlinear weak form F(u,p;v,q)
if n_glen == 1.0:
    print('setting-up weak form in unregularized Newtonian case (n_glen = 1.0) ...')
    F = ( inner(Bn * D(u), D(v)) - p * div(v) - div(u) * q \
          - inner(f_body, v) ) * dx \
        - inner(outflow_sigma, v) * ds(outflow_id)
else:
    print('using n_glen = %.3f and viscosity regularization eps = %.6f ...' \
          % (n_glen,args.eps))
    print('setting-up weak form ...')
    Du2 = 0.5 * inner(D(u), D(u)) + (args.eps * Dtyp)**2.0
    rr = 1.0/n_glen - 1.0
    F = ( inner(Bn * Du2**(rr/2.0) * D(u), D(v)) - p * div(v) - div(u) * q \
          - inner(f_body, v) ) * dx \
        - inner(outflow_sigma, v) * ds(outflow_id)

# slab-on-slope inflow boundary condition
C = (2.0 / (n_glen + 1.0)) * (rho * g * sin(args.alpha) / Bn)**n_glen
inflow_u = as_vector([C * (Hin**(n_glen+1.0) - (hsurfin - z)**(n_glen+1.0)), 0.0])

bcs = [ DirichletBC(Z.sub(0), noslip, base_id),
        DirichletBC(Z.sub(0), inflow_u, (inflow_id,)) ]

# solve
print('solving nonlinear variational problem ...')
solve(F == 0, up, bcs=bcs,
      options_prefix='s',
      solver_parameters={"snes_converged_reason": True,
                         #"ksp_converged_reason": True,
                         #"ksp_monitor": True,
                         "ksp_type": "fgmres",  # or "gmres" or "minres"
                         "pc_type": "fieldsplit",
                         "pc_fieldsplit_type": "schur",
                         "pc_fieldsplit_schur_factorization_type": "full",  # or "diag"
                         "fieldsplit_0_ksp_type": "preonly",
                         "fieldsplit_0_pc_type": "lu",
                         #"fieldsplit_0_ksp_converged_reason": True,
                         #"fieldsplit_1_ksp_converged_reason": True,
                         "fieldsplit_1_ksp_rtol": 1.0e-3,
                         "fieldsplit_1_ksp_type": "gmres",
                         "fieldsplit_1_pc_type": "none"})
# ALSO CONSIDER:
#    "ksp_type": "minres", "pc_type": "jacobi",
#    "mat_type": "aij", "ksp_type": "preonly", "pc_type": "svd",  # fully-direct solver
#    "mat_type": "aij", "ksp_view_mat": ":foo.m:ascii_matlab"

# report on and save solution parts
u,p = up.split()
u.rename('velocity')
p.rename('pressure')
# examine values directly:   print(p.vector().array())

P1 = FunctionSpace(mesh, "CG", 1)
one = Constant(1.0, domain=mesh)
area = assemble(dot(one,one) * dx)
#print('domain area = %.2e m2' % area)
pav = assemble(sqrt(dot(p, p)) * dx) / area
print('average pressure = %.2f Pa' % pav)
velmagav = assemble(sqrt(dot(u, u)) * dx) / area
print('average velocity magnitude = %.2f m a-1' % (secpera * velmagav))
umagmax = interpolate(sqrt(dot(u,u)),P1).dat.data.max()
print('maximum velocity magnitude = %.2f m a-1' % (secpera * umagmax))

# compute infinity-norm numerical errors relative to slab-on-slope when bs==0.0
if isslab:
    up_exact = Function(Z)
    u_exact,p_exact = up_exact.split()
    u_exact.interpolate(inflow_u)
    p_exact.interpolate(rho * g * cos(args.alpha) * (Hin - z))
    uerr = interpolate(sqrt(dot(u_exact-u,u_exact-u)),P1).dat.data.max()
    perr = interpolate(sqrt(dot(p_exact-p,p_exact-p)),W).dat.data.max()
    print('numerical errors: |u-uex|_inf = %.2e m a-1, |p-pex|_inf = %.2e Pa' \
          % (uerr*secpera,perr))

# write the solution to a file for visualisation with paraview
print('writing (velocity,pressure) to %s ...' % outname)
File(outname).write(u,p)

