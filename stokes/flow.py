#!/usr/bin/env python3
# (C) 2018--2020 Ed Bueler

# Solve glacier Glen-Stokes problem with explicit evolution of the free surface.
# See README.md for usage.

# TODO explore and package solver options

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from firedrake import *
from gendomain import Hin, L, bdryids, getdomaindims
from momentummodel import mixFEchoices, secpera, dayspera, MomentumModel
from meshmotion import surfacekinematical, movemesh
from surfaceutils import surfaceplot

# process options
parser = ArgumentParser(\
    description='Solve 2D glacier Glen-Stokes problem with evolving surface.  Requires an\nactivated Firedrake environment.  The Stokes solver has option prefix -s_...',
    formatter_class=RawTextHelpFormatter,add_help=False)
parser.add_argument('-alpha', type=float, default=0.1, metavar='X',
                    help='downward slope of bed as angle in radians (default=0.1)')
parser.add_argument('-deltat', type=float, default=0.0, metavar='X',
                    help='duration in days of time step for surface evolution\n(default=0.0)')
parser.add_argument('-Dtyp', type=float, default=2.0, metavar='X',
                    help='regularize viscosity using "+(eps Dtyp)^2" (default=2.0 a-1)')
                    # e.g. (800 m a-1) / 400 m = 2 a-1
parser.add_argument('-elements', metavar='X', default='P2P1',
                    choices=mixFEchoices,
                    help='stable mixed finite elements from: %s\n(default=P2P1)' \
                         % (','.join(mixFEchoices)) )
parser.add_argument('-eps', type=float, default=0.01, metavar='X',
                    help='regularize viscosity using "+(eps Dtyp)^2" (default=0.01)')
parser.add_argument('-flowhelp', action='store_true', default=False,
                    help='print help for flow.py options and stop')
parser.add_argument('-m', type=int, default=0, metavar='X',
                    help='number of time steps of deltat days of surface evolution\n(default=0)')
parser.add_argument('-mesh', metavar='MESH', type=str, default='',
                    help='input file name ending with .msh')
parser.add_argument('-n_glen', type=float, default=3.0, metavar='X',
                    help='Glen flow law exponent (default=3.0)')
parser.add_argument('-o', metavar='NAME', type=str, default='',
                    help='output file name ending with .pvd (default=MESH-.msh+.pvd)')
parser.add_argument('-osurface', metavar='NAME', type=str, default='',
                    help='save plot of surface values in this image file\n(i.e. .png,.pdf,...)')
parser.add_argument('-refine', type=int, default=0, metavar='N',
                    help='number of refinement levels after reading mesh (default=0)')
parser.add_argument('-save_rank', action='store_true',
                    help='add fields (element_rank,vertex_rank) to output file', default=False)
parser.add_argument('-sequence', type=int, default=0, metavar='N',
                    help='number of grid-sequencing levels (default=0)')
args, unknown = parser.parse_known_args()

if args.flowhelp:
    parser.print_help()
    sys.exit(0)
if len(args.mesh) == 0:
    print('ERROR: option -mesh required\n')
    parser.print_help()
    sys.exit(1)

# time-stepping requires both positive deltat and positive m, and does not
#     allow grid-sequencing
if args.m > 0:
    assert (args.deltat > 0.0)
if args.deltat > 0.0:
    assert (args.m > 0)
    assert (args.sequence == 0)

# output file name: strip .msh and replace with .pvd
if len(args.o) > 0:
    outname = args.o
else:
    outname = '.'.join(args.mesh.split('.')[:-1]) + '.pvd'

def printpar(thestr,comm=COMM_WORLD,indent=0):
    spaces = indent * '  '
    PETSc.Sys.Print('%s%s' % (spaces,thestr),comm=comm)

def describe(thismesh,indent=0):
    if thismesh.comm.size == 1:
        printpar('mesh has %d elements (cells) and %d vertices' \
                 % (thismesh.num_cells(),thismesh.num_vertices()),
                 indent=indent)
    else:
        PETSc.Sys.syncPrint('rank %d owns %d elements (cells) and can access %d vertices' \
                            % (thismesh.comm.rank,thismesh.num_cells(),thismesh.num_vertices()),
                            comm=thismesh.comm)
        PETSc.Sys.syncFlush(comm=thismesh.comm)

if args.sequence > 0:
    printpar('using %d levels of grid-sequencing ...' % args.sequence)

# read initial mesh, refine, and report on it
printpar('reading initial mesh from %s ...' % args.mesh,indent=args.sequence)
mesh = Mesh(args.mesh)
if args.refine + args.sequence > 0:
    hierarchy = MeshHierarchy(mesh, args.refine + args.sequence)
    mesh = hierarchy[args.refine]
describe(mesh,indent=args.sequence+1)
bs,bmin_initial,Hout_initial = getdomaindims(mesh)
printpar('geometry [m]: L = %.3f, bs = %.3f, Hin = %.3f' \
         %(L,bs,Hin),indent=args.sequence+1)
if bs < 1.0:
    printpar('slab geometry case ...',indent=args.sequence+2)
mesh._topology_dm.viewFromOptions('-dm_view')

# -osurface is not available in parallel
assert (mesh.comm.size == 1 or len(args.osurface) == 0)

# initialize momentum model
mm = MomentumModel()
mm.set_n_glen(args.n_glen)
mm.set_eps(args.eps)
mm.set_alpha(args.alpha)
mm.set_Dtyp_pera(args.Dtyp)
mm.set_Hin(Hin)
mm.set_Hout(Hout_initial)

outfile = File(outname)

# return integer-valued fields on the mesh which give the process rank for
# element ownership (piecewise constant; discontinuous) and vertex/node
# ownership (integer-valued on vertices; piecewise linear on mesh)
def getranks():
    element_rank = Function(FunctionSpace(mesh,'DG',0))
    element_rank.dat.data[:] = mesh.comm.rank
    element_rank.rename('element_rank')
    vertex_rank = Function(FunctionSpace(mesh,'CG',1))
    vertex_rank.dat.data[:] = mesh.comm.rank
    vertex_rank.rename('vertex_rank')
    return (element_rank,vertex_rank)

# solver mode: momentum-only solve of Stokes problem
def momentumsolve(ucoarse = None, pcoarse = None, indent=0):
    u,p = mm.solve(mesh,bdryids,args.elements,\
                   ucoarse=ucoarse,pcoarse=pcoarse)
    umagav,umagmax,pav,pmax = mm.solutionstats(mesh)
    printpar('flow speed: av = %10.3f m a-1,  max = %10.3f m a-1' \
             % (secpera*umagav,secpera*umagmax),
             indent=indent+1)
    printpar('pressure:   av = %10.3f bar,    max = %10.3f bar' \
             % (1.0e-5*pav,1.0e-5*pmax),
             indent=indent+1)
    return u,p

# solver mode: time-stepping loop with evolving surface
def timestepping():
    t_days = 0.0
    if args.save_rank:
        (element_rank,vertex_rank) = getranks()
    for j in range(args.m):
        printpar('solving at step %d: t = %.3f days ...' % (j,t_days))
        # get velocity, pressure, and mesh vertical displacement
        u,p = momentumsolve()
        dt = args.deltat * (secpera/dayspera)
        r = surfacekinematical(mesh,bdryids,u,dt)
        # save current state
        nu = mm.effectiveviscosity(mesh)
        if args.save_rank:
            outfile.write(u,p,nu,r,element_rank,vertex_rank, time=t_days)
        else:
            outfile.write(u,p,nu,r, time=t_days)
        # actually move mesh
        unstable = movemesh(mesh,r,bmin_initial)
        if unstable:
            printpar('\n\nSURFACE ELEVATION INSTABILITY DETECTED ... stopping\n\n')
            break
        # report on amount of movement
        t_days += args.deltat
        _,_,Hout = getdomaindims(mesh)
        mm.set_Hout(Hout)
        with r.dat.vec_ro as vr:
            absrmax = vr.norm(norm_type=PETSc.NormType.NORM_INFINITY)
        printpar('max. vert. disp. = %.3f m,  Hout = %.3f m' % (absrmax,Hout),
                 indent=1)
    return u,p,r,t_days

# in slab-on-slope case, compute numerical errors
def numericalerrorsslab(indent=0):
    if bs < 1.0:
        uerrmax,perrmax = mm.numerical_errors_slab(mesh)
        printpar('numerical errors: |u-uex|_inf = %.3e m a-1, |p-pex|_inf = %.3e Pa' \
                 % (uerrmax*secpera,perrmax),indent=indent)

# solve, after deciding on solver mode
if args.deltat > 0.0:
    printpar('writing (velocity,pressure,eff.visc.,vert.disp.) at each time step to %s ...' % outname)
    u,p,r,t_days = timestepping()  # returns at final time
elif args.sequence > 0:
    l = args.sequence
    printpar('solving for velocity and pressure ...',indent=l)
    u,p = momentumsolve(indent=l)
    for j in range(l):
        numericalerrorsslab(indent=l-j)
        mesh = hierarchy[args.refine+j+1]
        describe(mesh,indent=l-j-1)
        printpar('solving for velocity and pressure ...',indent=l-j-1)
        u,p = momentumsolve(ucoarse=u.copy(),pcoarse=p.copy(),indent=l-j-1)
else:
    printpar('solving for velocity and pressure ...')
    u,p = momentumsolve()
numericalerrorsslab()

# save results in paraview-readable file
nu = mm.effectiveviscosity(mesh)
if args.deltat > 0.0:
    printpar('writing END step (t = %.3f days) to file %s ...' % (t_days,outname))
    if args.save_rank:
        printpar('  writing (element_rank,vertex_rank) ...')
        outfile.write(u,p,nu,r,element_rank,vertex_rank, time=t_days)
    else:
        outfile.write(u,p,nu,r, time=t_days)
else:
    printpar('writing (velocity,pressure,eff.visc.) to %s ...' % outname)
    if args.save_rank:
        (element_rank,vertex_rank) = getranks()
        printpar('  writing (element_rank,vertex_rank) ...')
        outfile.write(u,p,nu,element_rank,vertex_rank)
    else:
        outfile.write(u,p,nu)

# generate image file with plot of surface values if desired
if len(args.osurface) > 0:
    if args.deltat > 0:
        surfaceplot(mesh,u,r,args.deltat,args.osurface)
    else:
        surfaceplot(mesh,u,None,args.deltat,args.osurface)

