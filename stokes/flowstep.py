#!/usr/bin/env python3
# (C) 2018 Ed Bueler

# Solve glacier bedrock-steps Glen-Stokes problem with evolving surface.  See
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
parser.add_argument('-deltat', type=float, default=0.0, metavar='X',
                    help='duration in days of time step for surface evolution')
parser.add_argument('-eps', type=float, default=0.01, metavar='X',
                    help='regularize viscosity using "+(eps Dtyp)^2" (default = 0.01)')
parser.add_argument('inname', metavar='INNAME',
                    help='input file name ending with .msh')
parser.add_argument('-m', type=int, default=1, metavar='X',
                    help='number of time steps of deltat days of surface evolution')
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
from genstepmesh import Hin, L, bdryids, getmeshdims
from physics import secpera, stokessolve, surfsolve, solutionstats, numericalerrorsslab

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

# extract mesh geometry needed in solver
bs,Hout = getmeshdims(mesh)
printpar('mesh geometry [m]: L = %.3f, bs = %.3f, Hin = %.3f, Hout = %.3f' \
         %(L,bs,Hin,Hout))
printpar('  bed slope angle alpha = %.6f radians' % args.alpha)
isslab = bs < 1.0
if isslab:
    printpar('  slab geometry case ...')

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

if args.n_glen == 1.0:
    printpar('unregularized Newtonian case (n_glen = 1.0)')
else:
    printpar('power law case with n_glen = %.3f and visc. reg. eps = %.6f' \
             % (args.n_glen,args.eps))

# time-stepping loop
t_days = 0.0
if args.deltat > 0.0:
    printpar('writing (velocity,pressure,displacement) at each time step to %s ...' % outname)
    outfile = File(outname)
up = Function(Z)
for j in range(args.m):
    if args.deltat > 0.0:
        printpar('step %d: t = %.3f days' % (j,t_days))
    printpar('  solving for velocity and pressure ...')

    # solve Stokes problem for u,p
    up = stokessolve(up,mesh,bdryids,Z,
                     Hin = Hin,
                     Hout = Hout,
                     n_glen = args.n_glen,
                     alpha = args.alpha,
                     eps = args.eps,
                     Dtyp = args.Dtyp / secpera)
    # report
    u,p = up.split()
    u.rename('velocity')
    p.rename('pressure')
    (umagav,umagmax,pav,pmax) = solutionstats(u,p,mesh)
    printpar('    flow speed: av = %10.3f m a-1,  max = %10.3f m a-1' \
             % (secpera*umagav,secpera*umagmax))
    printpar('    pressure:   av = %10.3f bar,    max = %10.3f bar' \
             % (1.0e-5*pav,1.0e-5*pmax))

    # time-stepping;  deltat=0 case is diagnostic only
    if args.deltat > 0.0:
        printpar('  solving for vertical mesh displacement from %.3f days motion ...' % args.deltat)
        r = surfsolve(mesh,bdryids,u)
        r *= args.deltat * (secpera/365.0)
        with r.dat.vec_ro as vr:
            absrmax = vr.norm(norm_type=PETSc.NormType.NORM_INFINITY)
        r.rename('displacement')
        printpar('    maximum vertical mesh displacement = %.3f m' % absrmax)
        printpar('  writing t = %.3f days values ...' % t_days)
        outfile.write(u,p,r, time=t_days)
        printpar('  applying displacement to mesh ...')
        Vc = mesh.coordinates.function_space()
        x,z = SpatialCoordinate(mesh)
        f = Function(Vc).interpolate(as_vector([x, z + r]))
        mesh.coordinates.assign(f)
        t_days += args.deltat

# compute numerical errors relative to slab-on-slope *if* bs==0.0
if isslab:
    (uerrmax,perrmax) = numericalerrorsslab(u,p,mesh,V,W,Hin,args.n_glen,args.alpha)
    printpar('numerical errors: |u-uex|_inf = %.3e m a-1, |p-pex|_inf = %.3e Pa' \
             % (uerrmax*secpera,perrmax))

if args.deltat > 0.0:
    # save final mesh
    printpar('step END: t = %.3f days' % t_days)
    printpar('  writing END values to finish %s ...' % outname)
    outfile.write(u,p,r, time=t_days)
else:
    # save solution if diagnostic
    printpar('writing (velocity,pressure) to %s ...' % outname)
    File(outname).write(u,p)

