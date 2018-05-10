#!/usr/bin/env python3
# (C) 2018 Ed Bueler

# Solve glacier bedrock-step Glen-Stokes problem.  See
# mccarthy/stokes/README.md for usage.

# process options
import argparse
mixFEchoices = ['P2P1','P3P2','P2P0','CRP0','P1P0']
parser = argparse.ArgumentParser(\
    description='Solve glacier Glen-Stokes problem on bedrock-step geometry.')
parser.add_argument('-alpha', type=float, default=0.1, metavar='X',
                    help='downward slope of bed as angle in radians (default = 0.1)')
parser.add_argument('-Dtyp', type=float, default=2.0, metavar='X',
                    help='regularize viscosity using "+(eps Dtyp)^2" (default = 2 a-1)')
                    # e.g. (800 m a-1) / 400 m = 2 a-1
parser.add_argument('-elements', metavar='X', default='P2P1',
                    choices=mixFEchoices,
                    help='stable mixed finite elements from: %s (default=P2P1)' \
                         % (','.join(mixFEchoices)) )
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

from firedrake import *
from firedrake.petsc import PETSc
from stepmesh import bdryids,getmeshdims
from physics import secpera,rho,g,getinflow,stokessolve

def printpar(thestr,comm=COMM_WORLD):
    PETSc.Sys.Print(thestr,comm=comm)

# read mesh and report on parallel decomposition (if appropriate)
printpar('reading initial mesh from %s ...' % args.inname)
mesh = Mesh(args.inname)
if mesh.comm.size == 1:
    printpar('  mesh has %d elements (cells) and %d vertices' \
          % (mesh.num_cells(),mesh.num_vertices()))
else:
    PETSc.Sys.syncPrint('  rank %d owns %d elements (cells) and can access %d vertices' \
                        % (mesh.comm.rank,mesh.num_cells(),mesh.num_vertices()), comm=mesh.comm)
    PETSc.Sys.syncFlush(comm=mesh.comm)

meshdims = getmeshdims(mesh)
printpar('mesh geometry [m]: L = %.3f, bs = %.3f, hsurfin = %.3f, Hout = %.3f' \
         %(meshdims['L'],meshdims['bs'],meshdims['hsurfin'],meshdims['Hout']))
printpar('using bed slope angle alpha = %.6f' % args.alpha)
if meshdims['isslab']:
    printpar('  ... detected slab geometry case ...')

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

# solve Stokes problem for u,p
printpar('setting-up weak form ...')
if args.n_glen == 1.0:
    printpar('  in unregularized Newtonian case (n_glen = 1.0) ...')
else:
    printpar('  using n_glen = %.3f and viscosity regularization eps = %.6f ...' \
          % (args.n_glen,args.eps))
printpar('solving for velocity and pressure ...')
up = stokessolve(mesh,bdryids,Z,meshdims,
                 n_glen=args.n_glen,
                 alpha=args.alpha,
                 eps=args.eps,
                 Dtyp=args.Dtyp / secpera)

# report on and save solution parts
u,p = up.split()
u.rename('velocity')
p.rename('pressure')

P1 = FunctionSpace(mesh, "CG", 1)
one = Constant(1.0, domain=mesh)
area = assemble(dot(one,one) * dx)
pav = assemble(sqrt(dot(p, p)) * dx) / area
printpar('average pressure = %.3f Pa' % pav)
velmagav = assemble(sqrt(dot(u, u)) * dx) / area
printpar('average velocity magnitude = %.3f m a-1' % (secpera * velmagav))
umag = interpolate(sqrt(dot(u,u)),P1)
with umag.dat.vec_ro as vumag:
    umagmax = vumag.max()[1]
printpar('maximum velocity magnitude = %.3f m a-1' % (secpera * umagmax))

# compute numerical errors relative to slab-on-slope *if* bs==0.0
if meshdims['isslab']:
    up_exact = Function(Z)
    u_exact,p_exact = up_exact.split()
    inflow_u = getinflow(mesh,meshdims,args.n_glen,args.alpha)
    u_exact.interpolate(inflow_u)
    Hin = meshdims['hsurfin'] - meshdims['bs']
    _,z = SpatialCoordinate(mesh)
    p_exact.interpolate(rho * g * cos(args.alpha) * (Hin - z))
    uerr = interpolate(sqrt(dot(u_exact-u,u_exact-u)),P1)
    perr = interpolate(sqrt(dot(p_exact-p,p_exact-p)),W)
    with uerr.dat.vec_ro as vuerr:
        uerrmax = vuerr.max()[1]
    with perr.dat.vec_ro as vperr:
        perrmax = vperr.max()[1]
    printpar('numerical errors: |u-uex|_inf = %.3e m a-1, |p-pex|_inf = %.3e Pa' \
          % (uerrmax*secpera,perrmax))

# write the solution to a file for visualisation with paraview
printpar('writing (velocity,pressure) to %s ...' % outname)
File(outname).write(u,p)

