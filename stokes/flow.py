#!/usr/bin/env python3
# (C) 2018--2020 Ed Bueler

# Solve glacier Glen-Stokes problem with evolving surface.  See README.md for usage.

# note this code misses out on geometric multigrid because it does not create a
# mesh hierarchy

# classic stability:
#   ./gendomain.py -bs 0.0 -o slab.geo
#   gmsh -2 slab.geo
#   ./flow.py -deltat 20.0 -m 50 slab.msh  # have run to -m 500 w evident stability
# classic surface instability; sawtooth emerges at center of top of mesh:
#   ./flow.py -deltat 40.0 -m 25 slab.msh
# apparent time step limit: about dt=35 days
#   (but instability flows off end; dt=37 causes full-thickness blowup)

# for grids refined by factor of 2 it matters *A LOT* how that refinement happens
#   use of "./gendomain.py -bs 0.0 -refine 2 -o slab2.geo"  gives dt=1 or worse
#   use of gmsh to split cells ONCE gives dt=15 or better

# process options
from argparse import ArgumentParser, RawTextHelpFormatter
mixFEchoices = ['P2P1','P3P2','P2P0','CRP0','P1P0']
parser = ArgumentParser(\
    description='Solve 2D glacier Glen-Stokes problem with evolving surface.  Requires an activated Firedrake environment.',
    formatter_class=RawTextHelpFormatter,add_help=False)
parser.add_argument('-flowhelp', action='store_true', default=False,
                    help='print help for flow.py options and stop')
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
parser.add_argument('-m', type=int, default=1, metavar='X',
                    help='number of time steps of deltat days of surface evolution')
parser.add_argument('-mesh', metavar='MESH', type=str, default='',
                    help='input file name ending with .msh')
parser.add_argument('-n_glen', type=float, default=3.0, metavar='X',
                    help='Glen flow law exponent (default = 3.0)')
parser.add_argument('-o', metavar='NAME', type=str, default='',
                    help='output file name ending with .pvd (default = MESH-.msh+.pvd)')
parser.add_argument('-save_rank', action='store_true',
                    help='add fields (element_rank,vertex_rank) to output file', default=False)
parser.add_argument('-osurface', metavar='NAME', type=str, default='',
                    help='save a plot of surface values of (h,u,w) in this image file (.png,.pdf,...)')
args, unknown = parser.parse_known_args()

import sys
if args.flowhelp:
    parser.print_help()
    sys.exit(0)
if len(args.mesh) == 0:
    print('ERROR: option -mesh required')
    print('')
    parser.print_help()
    sys.exit(0)

if len(args.o) > 0:
    outname = args.o
else:
    outname = '.'.join(args.mesh.split('.')[:-1]) + '.pvd'  # strip .msh and replace with .pvd

from firedrake import *
from firedrake.petsc import PETSc
from gendomain import Hin, L, bdryids, getdomaindims
from physics import secpera, stokessolve, solutionstats, numericalerrors_slab
from meshactions import getranks, getsurfaceelevation, solvevdisp, \
                        getsurfacevdispfunction, getsurfacevelocityfunction

def printpar(thestr,comm=COMM_WORLD):
    PETSc.Sys.Print(thestr,comm=comm)

# read mesh and report on parallel decomposition (if appropriate)
printpar('reading initial mesh from %s ...' % args.mesh)
mesh = Mesh(args.mesh)
assert (mesh.comm.size == 1 or len(args.osurface) == 0)  # -osurface not available in parallel
if mesh.comm.size == 1:
    printpar('  mesh has %d elements (cells) and %d vertices' \
          % (mesh.num_cells(),mesh.num_vertices()))
else:
    PETSc.Sys.syncPrint('  rank %d owns %d elements (cells) and can access %d vertices' \
                        % (mesh.comm.rank,mesh.num_cells(),mesh.num_vertices()), comm=mesh.comm)
    PETSc.Sys.syncFlush(comm=mesh.comm)

# extract mesh geometry needed in solver
bs,bmin_initial,Hout = getdomaindims(mesh)
printpar('geometry [m]: L = %.3f, bs = %.3f, Hin = %.3f, Hout = %.3f' \
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
    printpar('writing (velocity,pressure,vertical_displacement) at each time step to %s ...' % outname)
    outfile = File(outname)
up = Function(Z)
P1 = FunctionSpace(mesh,'CG',1)  # used for mesh displacement
if args.save_rank:
    (element_rank,vertex_rank) = getranks(mesh)
for j in range(args.m):
    if args.deltat > 0.0:
        printpar('step %d: t = %.3f days' % (j,t_days))
    printpar('  solving for velocity and pressure ...')

    # solve Stokes problem for u,p; defined in physics.py; prefix s_
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
        printpar('  solving kinematical equation for vertical mesh displacement rate ...')
        # use surface kinematical equation to get boundary condition for mesh displacement problem
        h = getsurfaceelevation(mesh,bdryids['top'])
        x,z = SpatialCoordinate(mesh)
        xval = Function(P1).interpolate(x)
        zval = Function(P1).interpolate(z)
        phi = Function(P1)
        phi.dat.data[:] = zval.dat.data_ro - h(xval.dat.data_ro)
        dt = args.deltat * (secpera/365.0)
        # FIXME add in climatic mass balance a(x) here; want h_t = a - u[0] h_x + u[1]
        #       currently uses:  a = Constant(0.0)
        deltah = Function(P1).interpolate( dt * (Constant(0.0) + dot(grad(phi),u)) )
        # solve mesh displacement problem; defined in meshactions.py; prefix vd_
        r = solvevdisp(mesh,bdryids,deltah)
        with r.dat.vec_ro as vr:
            absrmax = vr.norm(norm_type=PETSc.NormType.NORM_INFINITY)
        r.rename('vertical_displacement')
        # save complete current state
        if args.save_rank:
            outfile.write(u,p,r,element_rank,vertex_rank, time=t_days)
        else:
            outfile.write(u,p,r, time=t_days)
        # actually move mesh
        Vc = mesh.coordinates.function_space()
        f = Function(Vc).interpolate(as_vector([x, z + r]))
        mesh.coordinates.assign(f)
        # report on amount of movement
        t_days += args.deltat
        bs,_,Hout = getdomaindims(mesh)
        printpar('    max. vert. mesh disp. = %.3f m,  Hout = %.3f m' % (absrmax,Hout))
        if any(f.dat.data_ro[:,1] < bmin_initial - 1.0):  # mesh z values below bed is extreme instability
            printpar('\n\nSURFACE ELEVATION INSTABILITY DETECTED ... stopping\n\n')
            break

# compute numerical errors relative to slab-on-slope *if* bs==0.0
if isslab:
    (uerrmax,perrmax) = numericalerrors_slab(u,p,mesh,V,W,Hin,args.n_glen,args.alpha)
    printpar('numerical errors: |u-uex|_inf = %.3e m a-1, |p-pex|_inf = %.3e Pa' \
             % (uerrmax*secpera,perrmax))

# save results in .pvd
if args.deltat > 0.0:
    # save final mesh
    printpar('step END: t = %.3f days' % t_days)
    printpar('  writing values to finish file %s ...' % outname)
    if args.save_rank:
        outfile.write(u,p,r,element_rank,vertex_rank, time=t_days)
    else:
        outfile.write(u,p,r, time=t_days)
else:
    # save solution if diagnostic
    printpar('writing (velocity,pressure) to %s ...' % outname)
    if args.save_rank:
        File(outname).write(u,p,element_rank,vertex_rank)
    else:
        File(outname).write(u,p)

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
    hfcn = getsurfaceelevation(mesh,bdryids['top'])
    ufcn,wfcn = getsurfacevelocityfunction(mesh,bdryids['top'],Z,u)
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
    plt.plot(x,secpera*ufcn(x),label='horizontal velocity')
    plt.ylabel('u  [m/a]')
    plt.legend()
    removexticks()
    plt.subplot(rows,1,3)
    plt.plot(x,secpera*wfcn(x),label='vertical velocity')
    plt.ylabel('w  [m/a]')
    plt.legend()
    if args.deltat > 0.0:
        removexticks()
        rfcn = getsurfacevdispfunction(mesh,bdryids['top'],r)
        plt.subplot(rows,1,4)
        plt.plot(x,rfcn(x)/(args.deltat/365.2422),'r',label='surface elevation rate (last time step)')
        plt.ylabel('h_t  [m/a]')
        plt.legend()
    plt.xlabel('x  [m]')
    plt.savefig(args.osurface,bbox_inches='tight')

