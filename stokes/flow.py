#!/usr/bin/env python3
# (C) 2018--2020 Ed Bueler

# Solve glacier Glen-Stokes problem with explicit evolution of the free surface.
# See README.md for usage.

# FIXME this code misses out on grid-sequencing and geometric multigrid because
#       it does not create a mesh hierarchy

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from firedrake import *
from gendomain import Hin, L, bdryids, getdomaindims
from momentummodel import mixFEchoices, secpera, dayspera, MomentumModel
from meshmotion import surfacekinematical, movemesh
from surfaceutils import getsurfaceelevation, surfaceplot

# process options
parser = ArgumentParser(\
    description='Solve 2D glacier Glen-Stokes problem with evolving surface.  Requires an activated Firedrake environment.',
    formatter_class=RawTextHelpFormatter,add_help=False)
parser.add_argument('-alpha', type=float, default=0.1, metavar='X',
                    help='downward slope of bed as angle in radians (default = 0.1)')
parser.add_argument('-deltat', type=float, default=0.0, metavar='X',
                    help='duration in days of time step for surface evolution')
parser.add_argument('-Dtyp', type=float, default=2.0, metavar='X',
                    help='regularize viscosity using "+(eps Dtyp)^2" (default = 2 a-1)')
                    # e.g. (800 m a-1) / 400 m = 2 a-1
parser.add_argument('-elements', metavar='X', default='P2P1',
                    choices=mixFEchoices,
                    help='stable mixed finite elements from: %s (default=P2P1)' \
                         % (','.join(mixFEchoices)) )
parser.add_argument('-eps', type=float, default=0.01, metavar='X',
                    help='regularize viscosity using "+(eps Dtyp)^2" (default = 0.01)')
parser.add_argument('-flowhelp', action='store_true', default=False,
                    help='print help for flow.py options and stop')
parser.add_argument('-m', type=int, default=1, metavar='X',
                    help='number of time steps of deltat days of surface evolution')
parser.add_argument('-mesh', metavar='MESH', type=str, default='',
                    help='input file name ending with .msh')
parser.add_argument('-n_glen', type=float, default=3.0, metavar='X',
                    help='Glen flow law exponent (default = 3.0)')
parser.add_argument('-o', metavar='NAME', type=str, default='',
                    help='output file name ending with .pvd (default = MESH-.msh+.pvd)')
parser.add_argument('-osurface', metavar='NAME', type=str, default='',
                    help='save a plot of surface values of (h,u,w) in this image file (.png,.pdf,...)')
parser.add_argument('-save_rank', action='store_true',
                    help='add fields (element_rank,vertex_rank) to output file', default=False)
args, unknown = parser.parse_known_args()

if args.flowhelp:
    parser.print_help()
    sys.exit(0)
if len(args.mesh) == 0:
    print('ERROR: option -mesh required\n')
    parser.print_help()
    sys.exit(1)

# output file name: strip .msh and replace with .pvd
if len(args.o) > 0:
    outname = args.o
else:
    outname = '.'.join(args.mesh.split('.')[:-1]) + '.pvd'

def printpar(thestr,comm=COMM_WORLD):
    PETSc.Sys.Print(thestr,comm=comm)

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

# read initial mesh and report on it
printpar('reading initial mesh from %s ...' % args.mesh)
mesh = Mesh(args.mesh)
if mesh.comm.size == 1:
    printpar('  mesh has %d elements (cells) and %d vertices' \
             % (mesh.num_cells(),mesh.num_vertices()))
else:
    PETSc.Sys.syncPrint('  rank %d owns %d elements (cells) and can access %d vertices' \
                        % (mesh.comm.rank,mesh.num_cells(),mesh.num_vertices()),
                        comm=mesh.comm)
    PETSc.Sys.syncFlush(comm=mesh.comm)
bs,bmin_initial,Hout_initial = getdomaindims(mesh)
printpar('  geometry [m]: L = %.3f, bs = %.3f, Hin = %.3f' \
         %(L,bs,Hin))
if bs < 1.0:
    printpar('    slab geometry case ...')
mesh._topology_dm.viewFromOptions('-dm_view')

# -osurface is not available in parallel
assert (mesh.comm.size == 1 or len(args.osurface) == 0)

# initialize momentum model, which will own the mesh from now on
mm = MomentumModel()
mm.set_n_glen(args.n_glen)
mm.set_eps(args.eps)
mm.set_alpha(args.alpha)
mm.set_Dtyp_pera(args.Dtyp)
mm.set_Hin(Hin)
mm.set_Hout(Hout_initial)

# time-stepping loop
outfile = File(outname)
if args.deltat > 0.0:
    printpar('writing (velocity,pressure,vertical_displacement) at each time step to %s ...' % outname)
if args.save_rank:
    printpar('writing (element_rank,vertex_rank) into output file %s ...' % outname)
    (element_rank,vertex_rank) = getranks(mesh)
t_days = 0.0
for j in range(args.m):
    if args.deltat > 0.0:
        printpar('step %d: t = %.3f days' % (j,t_days))

    # solve Stokes problem; solver prefix s_
    printpar('solving for velocity and pressure ...')
    u,p = mm.solve(mesh,bdryids,args.elements)

    # report
    umagav,umagmax,pav,pmax = mm.solutionstats(mesh)
    printpar('  flow speed: av = %10.3f m a-1,  max = %10.3f m a-1' \
             % (secpera*umagav,secpera*umagmax))
    printpar('  pressure:   av = %10.3f bar,    max = %10.3f bar' \
             % (1.0e-5*pav,1.0e-5*pmax))

    # time-stepping;  deltat=0 case is diagnostic only
    if args.deltat > 0.0:
        printpar('solving kinematical equation for vertical mesh displacement rate ...')
        dt = args.deltat * (secpera/dayspera)
        # r = mesh vertical displacement
        r = surfacekinematical(mesh,bdryids,u,dt)
        # save complete current state
        if args.save_rank:
            outfile.write(u,p,r,element_rank,vertex_rank, time=t_days)
        else:
            outfile.write(u,p,r, time=t_days)
        # actually move mesh
        mesh,unstable = movemesh(mesh,r,bmin_initial)
        if unstable:
            printpar('\n\nSURFACE ELEVATION INSTABILITY DETECTED ... stopping\n\n')
            break
        # report on amount of movement
        t_days += args.deltat
        _,_,Hout = getdomaindims(mesh)
        mm.set_Hout(Hout)
        with r.dat.vec_ro as vr:
            absrmax = vr.norm(norm_type=PETSc.NormType.NORM_INFINITY)
        printpar('  max. vert. mesh disp. = %.3f m,  Hout = %.3f m' % (absrmax,Hout))

# compute numerical errors relative to slab-on-slope *if* bs==0.0
if bs < 1.0:
    uerrmax,perrmax = mm.numerical_errors_slab(mesh)
    printpar('numerical errors: |u-uex|_inf = %.3e m a-1, |p-pex|_inf = %.3e Pa' \
             % (uerrmax*secpera,perrmax))

# save results in .pvd
if args.deltat > 0.0:
    printpar('writing END step (t = %.3f days) to file %s ...' % (t_days,outname))
    if args.save_rank:
        outfile.write(u,p,r,element_rank,vertex_rank, time=t_days)
    else:
        outfile.write(u,p,r, time=t_days)
else:
    printpar('writing (velocity,pressure) to %s ...' % outname)
    if args.save_rank:
        outfile.write(u,p,element_rank,vertex_rank)
    else:
        outfile.write(u,p)

# generate plot of surface values if desired
if len(args.osurface) > 0:
    surfaceplot(mesh,u,r,args.deltat,args.osurface)

