#!/usr/bin/env python3
# (C) 2018--2020 Ed Bueler

# Solve glacier Glen-Stokes problem with explicit evolution of the free surface.
# See README.md for usage.

# FIXME this code misses out on geometric multigrid because it does not create a
# mesh hierarchy

# process options
from argparse import ArgumentParser, RawTextHelpFormatter
from momentummodel import mixFEchoices, MomentumModel

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

import sys
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

from firedrake import *
from gendomain import Hin, L, bdryids, getdomaindims
from meshactions import getranks, getsurfaceelevation, solvevdisp, \
                        getsurfacevdispfunction, getsurfacevelocityfunction

def printpar(thestr,comm=COMM_WORLD):
    PETSc.Sys.Print(thestr,comm=comm)

# read initial mesh and report on parallel decomposition (if appropriate)
printpar('reading initial mesh from %s ...' % args.mesh)
mesh0 = Mesh(args.mesh)
assert (mesh0.comm.size == 1 or len(args.osurface) == 0)  # -osurface not available in parallel
if mesh0.comm.size == 1:
    printpar('  mesh has %d elements (cells) and %d vertices' \
          % (mesh0.num_cells(),mesh0.num_vertices()))
else:
    PETSc.Sys.syncPrint('  rank %d owns %d elements (cells) and can access %d vertices' \
                        % (mesh0.comm.rank,mesh0.num_cells(),mesh0.num_vertices()), comm=mesh0.comm)
    PETSc.Sys.syncFlush(comm=mesh0.comm)

# extract mesh geometry needed in solver
bs,bmin_initial,Hout_initial = getdomaindims(mesh0)
printpar('  geometry [m]: L = %.3f, bs = %.3f, Hin = %.3f' \
         %(L,bs,Hin))
isslab = bs < 1.0
if isslab:
    printpar('    slab geometry case ...')

# initialize momentum model
mm = MomentumModel(mesh0,bdryids,args.elements)
mm.set_nglen(args.n_glen)
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
    (element_rank,vertex_rank) = getranks(mm.mesh)  # FIXME mm.mesh should be private
t_days = 0.0
for j in range(args.m):
    if args.deltat > 0.0:
        printpar('step %d: t = %.3f days' % (j,t_days))

    # solve Stokes problem; defined in FIXME physics.py; prefix s_
    printpar('solving for velocity and pressure ...')
    mm.solve()

    # report
    umagav,umagmax,pav,pmax = mm.solutionstats()
    printpar('  flow speed: av = %10.3f m a-1,  max = %10.3f m a-1' \
             % (mm.secpera()*umagav,mm.secpera()*umagmax))
    printpar('  pressure:   av = %10.3f bar,    max = %10.3f bar' \
             % (1.0e-5*pav,1.0e-5*pmax))

    # time-stepping;  deltat=0 case is diagnostic only
    if args.deltat > 0.0:
        printpar('solving kinematical equation for vertical mesh displacement rate ...')
        # use surface kinematical equation to get boundary condition for mesh displacement problem
        # FIXME mm.mesh should be private
        h = getsurfaceelevation(mm.mesh,bdryids['top'])
        x,z = SpatialCoordinate(mm.mesh)
        P1 = FunctionSpace(mm.mesh,'CG',1)
        xval = Function(P1).interpolate(x)
        zval = Function(P1).interpolate(z)
        phi = Function(P1)
        phi.dat.data[:] = zval.dat.data_ro - h(xval.dat.data_ro)
        dt = args.deltat * (mm.secpera()/365.0)  # FIXME inconsistent with below
        # FIXME add in climatic mass balance a(x) here; want h_t = a - u[0] h_x + u[1]
        #       currently uses:  a = Constant(0.0)
        u,p = mm.getsolution()
        deltah = Function(P1).interpolate( dt * (Constant(0.0) + dot(grad(phi),u)) )
        # solve mesh displacement problem; defined in meshactions.py; prefix vd_
        r = solvevdisp(mm.mesh,bdryids,deltah)
        with r.dat.vec_ro as vr:
            absrmax = vr.norm(norm_type=PETSc.NormType.NORM_INFINITY)
        r.rename('vertical_displacement')
        # save complete current state
        if args.save_rank:
            outfile.write(u,p,r,element_rank,vertex_rank, time=t_days)
        else:
            outfile.write(u,p,r, time=t_days)
        # actually move mesh
        Vc = mm.mesh.coordinates.function_space()
        f = Function(Vc).interpolate(as_vector([x, z + r]))
        mm.set_mesh_coordinates(f)
        # report on amount of movement
        t_days += args.deltat
        bs,_,Hout = getdomaindims(mm.mesh)
        mm.set_Hout(Hout)
        printpar('  max. vert. mesh disp. = %.3f m,  Hout = %.3f m' % (absrmax,Hout))
        if any(f.dat.data_ro[:,1] < bmin_initial - 1.0):  # mesh z values below bed is extreme instability
            printpar('\n\nSURFACE ELEVATION INSTABILITY DETECTED ... stopping\n\n')
            break

# compute numerical errors relative to slab-on-slope *if* bs==0.0
if isslab:
    uerrmax,perrmax = mm.numerical_errors_slab()
    printpar('numerical errors: |u-uex|_inf = %.3e m a-1, |p-pex|_inf = %.3e Pa' \
             % (uerrmax*mm.secpera(),perrmax))

# save results in .pvd
u,p = mm.getsolution()  # FIXME?
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

def removexticks():
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

# generate plot of surface values if desired
if len(args.osurface) > 0:
    import numpy as np
    import matplotlib.pyplot as plt
    x = np.linspace(0.0,L,401)
    hfcn = getsurfaceelevation(mm.mesh,bdryids['top'])
    ufcn,wfcn = getsurfacevelocityfunction(mm.mesh,bdryids['top'],mm.Z,u)
    plt.figure(figsize=(6.0,8.0))
    if args.deltat > 0.0:
        rows = 4
        printpar('plotting surface values of (h,u,w,h_t) in file %s ...' % args.osurface)
    else:
        rows = 3
        printpar('plotting surface values of (h,u,w) in file %s ...' % args.osurface)
    plt.subplot(rows,1,1)
    plt.plot(x,hfcn(x),'g',label='surface elevation')
    plt.ylabel('h  [m]')
    plt.legend()
    removexticks()
    plt.subplot(rows,1,2)
    plt.plot(x,mm.secpera()*ufcn(x),label='horizontal velocity')
    plt.ylabel('u  [m/a]')
    plt.legend()
    removexticks()
    plt.subplot(rows,1,3)
    plt.plot(x,mm.secpera()*wfcn(x),label='vertical velocity')
    plt.ylabel('w  [m/a]')
    plt.legend()
    if args.deltat > 0.0:
        removexticks()
        rfcn = getsurfacevdispfunction(mm.mesh,bdryids['top'],r)
        plt.subplot(rows,1,4)
        plt.plot(x,rfcn(x)/(args.deltat/365.2422),'r',label='surface elevation rate (last time step)')
        plt.ylabel('h_t  [m/a]')
        plt.legend()
    plt.xlabel('x  [m]')
    plt.savefig(args.osurface,bbox_inches='tight')

