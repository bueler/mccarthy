#!/usr/bin/env python3
# (C) 2018--2024 Ed Bueler

# Solves glacier Glen-Stokes problem.  See README.md for usage.
# Documented by doc/doc.pdf.

# process options
packagechoices = ['SchurDirect','Direct']
from argparse import ArgumentParser, RawTextHelpFormatter
parser = ArgumentParser(description=\
'''Solve 2D glacier Glen-Stokes problem.  Requires an activated
Firedrake environment.  Solver uses option prefix -s_...''',
    formatter_class=RawTextHelpFormatter, add_help=False)
parser.add_argument('-alpha', type=float, default=0.1, metavar='X',
                    help='downward slope of bed as angle in radians (default=0.1)')
parser.add_argument('-Dtyp', type=float, default=2.0, metavar='X',
                    help='regularize viscosity using "+(eps Dtyp)^2" (default=2.0 a-1)')
                    # e.g. (800 m a-1) / 400 m = 2 a-1
parser.add_argument('-eps', type=float, default=0.01, metavar='X',
                    help='regularize viscosity using "+(eps Dtyp)^2" (default=0.01)')
parser.add_argument('-flowhelp', action='store_true', default=False,
                    help='print help for flow.py options and exit')
parser.add_argument('-Hin', type=float, default=400.0, metavar='X',
                    help='upstream thickness of ice (default=400 m)')
parser.add_argument('-Hout', type=float, default=400.0, metavar='X',
                    help='downstream thickness of ice (default=400 m)')
parser.add_argument('-mesh', metavar='MESH', type=str, default='',
                    help='input file name ending with .msh')
parser.add_argument('-o', metavar='NAME', type=str, default='',
                    help='output file name ending with .pvd')
parser.add_argument('-osurface', metavar='NAME', type=str, default='',
                    help='save plot of surface values in this image file\n(i.e. .png,.pdf,...)')
parser.add_argument('-package', metavar='X', default='Direct',
                    choices=packagechoices,
                    help='solver package from:\n%s\n(default=Direct)' \
                         % (','.join(packagechoices)) )
parser.add_argument('-refine', type=int, default=0, metavar='N',
                    help='number of refinement levels after reading mesh (default=0)')
parser.add_argument('-sequence', type=int, default=0, metavar='N',
                    help='number of grid-sequencing levels (default=0)')
parser.add_argument('-slab', action='store_true',
                    help='compute errors relative to slab-on-slope')
args, passthroughoptions = parser.parse_known_args()

import sys
if args.flowhelp:
    parser.print_help()
    sys.exit(0)

import petsc4py
petsc4py.init(passthroughoptions)
from firedrake import *
from firedrake.output import VTKFile
from firedrake.petsc import PETSc

import numpy as np
from domain import bdryids
from momentummodel import secpera, dayspera, MomentumModel
from surfaceutils import surfaceplot

if len(args.mesh) == 0:
    print('ERROR: option -mesh required\n')
    parser.print_help()
    sys.exit(1)

def printpar(thestr, comm=COMM_WORLD, indent=0):
    spaces = indent * '  '
    PETSc.Sys.Print('%s%s' % (spaces, thestr), comm=comm)

def describe(thismesh, indent=0):
    if thismesh.comm.size == 1:
        printpar('mesh has %d elements (cells) and %d vertices' \
                 % (thismesh.num_cells(),thismesh.num_vertices()),
                 indent=indent)
    else:
        PETSc.Sys.syncPrint('rank %d owns %d elements (cells) and can access %d vertices' \
                            % (thismesh.comm.rank,thismesh.num_cells(),thismesh.num_vertices()), comm=thismesh.comm)
        PETSc.Sys.syncFlush(comm=thismesh.comm)

if args.sequence > 0:
    printpar('using %d levels of grid-sequencing ...' % args.sequence)

# read initial mesh, refine, and report on it
printpar('reading initial mesh from %s ...' % args.mesh,indent=args.sequence)
mesh = Mesh(args.mesh)
if args.refine + args.sequence > 0:
    hierarchy = MeshHierarchy(mesh, args.refine + args.sequence)
    if args.refine > 0:
        printpar('refining mesh %d times ...' % args.refine,indent=args.sequence+1)
    mesh = hierarchy[args.refine]
describe(mesh, indent=args.sequence+1)

mesh.topology_dm.viewFromOptions('-dm_view')

# -osurface is not available in parallel
assert (mesh.comm.size == 1 or len(args.osurface) == 0)

# initialize momentum model
mm = MomentumModel(eps=args.eps, alpha=args.alpha,
                   Dtyp_pera=args.Dtyp, Hin=args.Hin, Hout=args.Hout)

# solver mode: momentum-only solve of Stokes problem
def momentumsolve(package=args.package, upold = None, upcoarse = None, indent=0):
    up = mm.solve(mesh, bdryids,
                  package=package, upold=upold, upcoarse=upcoarse)
    umagav,umagmax,pav,pmax = mm.solutionstats(mesh)
    printpar('flow speed: av = %10.3f m a-1,  max = %10.3f m a-1' \
             % (secpera*umagav,secpera*umagmax),
             indent=indent+1)
    printpar('pressure:   av = %10.3f bar,    max = %10.3f bar' \
             % (1.0e-5*pav,1.0e-5*pmax),
             indent=indent+1)
    return up

# in slab-on-slope case, compute and report numerical errors
def numericalerrorsslab(indent=0):
    uexactnorm, uerr, pexactnorm, perr = mm.numerical_errors_slab(mesh)
    printpar('numerical errors: |u-uex|_2/|uex|_2 = %.3e, |p-pex|_2/|pex|_2 = %.3e' \
             % (uerr / uexactnorm, perr / pexactnorm), indent=indent)

# solve, after deciding on solver mode
printpar('using %s solver package ...' % args.package)
if args.sequence > 0:
    l = args.sequence
    printpar(f'solving for velocity and pressure using {l} levels of grid-sequencing ...', indent=l)
    up = momentumsolve(package=args.package, indent=l)
    for j in range(l):
        if args.slab:
            numericalerrorsslab(indent=l-j)
        mesh = hierarchy[args.refine+j+1]
        printpar('transferring to next-refined mesh ...', indent=l-j-1)
        describe(mesh, indent=l-j-1)
        printpar('solving for velocity and pressure ...', indent=l-j-1)
        up = momentumsolve(package=args.package, upcoarse=up, indent=l-j-1)
else:
    printpar('solving for velocity and pressure ...')
    up = momentumsolve(package=args.package)
if args.slab:
    numericalerrorsslab()
u, p = up.subfunctions

# save results in paraview-readable file
if len(args.o) > 0:
    nu = mm.effectiveviscosity(mesh)
    if mesh.comm.size > 1:
        rank = Function(FunctionSpace(mesh,'DG',0))
        rank.dat.data[:] = mesh.comm.rank
        rank.rename('rank (element-wise)')
        printpar(f'writing (velocity, pressure, eff.visc., rank) to {args.o} ...')
        VTKFile(args.o).write(u, p, nu, rank)
    else:
        printpar(f'writing (velocity,pressure,eff.visc.) to {args.o} ...')
        VTKFile(args.o).write(u, p, nu)

# generate image file with plot of surface values if desired
if len(args.osurface) > 0:
    surfaceplot(mesh, u, args.osurface)
