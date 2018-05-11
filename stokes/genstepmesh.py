#!/usr/bin/env python3
# (C) 2018 Ed Bueler

# Usage:
#   $ ./genstepmesh.py glacier.geo
#   $ gmsh -2 glacier.geo
#   ...
# which generates glacier.msh.  (Use "gmsh glacier.geo" for GUI version.)
# One may inspect the mesh by loading it in python/firedrake:
#   $ source firedrake/bin/activate
#   (firedrake) $ ipython
#   ...
#   In [1]: from firedrake import *
#   In [2]: Mesh('glacier.msh')
# Main purpose is to solve the Stokes equations on this domain.  See
# README.md and flowstep.py.

# immutable domain distances in meters
Hin = 400.0      # input (and initial output) thickness (z)
L = 2000.0       # total along-flow length (x)
Ls = 1300.0      # location of bedrock step (x)
Ks = 150.0       # distance over which initial transition in surface happens (x)

# numbering of parts of boundary *must match generation script genstepmesh.py*
bdryids = {'outflow' : 41,
           'top'     : 42,
           'inflow'  : 43,
           'base'    : 44}

def writegeometry(geo,bs):
    # points on boundary
    Hout = Hin
    hsurfin = Hin + bs
    geo.write('Point(1) = {%f,%f,0,lc};\n' % (L,0.0))
    geo.write('Point(2) = {%f,%f,0,lc};\n' % (L,Hout))
    geo.write('Point(3) = {%f,%f,0,lc};\n' % (Ls+Ks,Hout))
    # FIXME  between Point(3) and Point(4) could be a nice smooth curve
    geo.write('Point(4) = {%f,%f,0,lc};\n' % (Ls-4*Ks,hsurfin))
    geo.write('Point(5) = {%f,%f,0,lc};\n' % (0.0,hsurfin))
    geo.write('Point(6) = {%f,%f,0,lc};\n' % (0.0,bs))
    if abs(bs) >= 1.0:
        geo.write('Point(7) = {%f,%f,0,lc_corner};\n' % (Ls,bs))
        geo.write('Point(8) = {%f,%f,0,lc};\n' % (Ls,0.0))

    # lines along boundary
    geo.write('Line(11) = {1,2};\n')
    geo.write('Line(12) = {2,3};\n')
    geo.write('Line(13) = {3,4};\n')
    geo.write('Line(14) = {4,5};\n')
    geo.write('Line(15) = {5,6};\n')
    if abs(bs) >= 1.0:
        geo.write('Line(16) = {6,7};\n')
        geo.write('Line(17) = {7,8};\n')
        geo.write('Line(18) = {8,1};\n')
    else:
        geo.write('Line(16) = {6,1};\n')

    # loop and surface allows defining a 2D mesh
    if abs(bs) >= 1.0:
        geo.write('Line Loop(21) = {11,12,13,14,15,16,17,18};\n')
    else:
        geo.write('Line Loop(21) = {11,12,13,14,15,16};\n')
    geo.write('Plane Surface(31) = {21};\n')

    # "Physical" for marking boundary conditions
    geo.write('Physical Line(%d) = {11};\n' % bdryids['outflow'])
    geo.write('Physical Line(%d) = {12,13,14};\n' % bdryids['top'])
    geo.write('Physical Line(%d) = {15};\n' % bdryids['inflow'])
    if abs(bs) >= 1.0:
        geo.write('Physical Line(%d) = {16,17,18};\n' % bdryids['base'])
    else:
        geo.write('Physical Line(%d) = {16};\n' % bdryids['base'])
    geo.write('Physical Surface(51) = {31};\n')   # ensure all interior elements written(?)

# dynamically extract mesh geometry making these definitions (tolerance=1cm):
#   bs      = height of bed step         = (min z-coordinate at x=0)
#   Hout    = ice thickness at output    = (max z-coordinate at x=L)
#   isslab  = (bs > 0 and Hout == Hin)
# in parallel no process owns the whole mesh so MPI_Allreduce() is needed
def getmeshdims(mesh,tol=0.01):
    from mpi4py import MPI
    dims = {}
    # mesh needs to be a Mesh from Firedrake
    xa = mesh.coordinates.dat.data_ro[:,0]  # .data_ro acts like VecGetArrayRead
    za = mesh.coordinates.dat.data_ro[:,1]
    loc_bs = 9.9999e99
    if any(xa < tol):              # some processes may have no such points
        loc_bs = min(za[xa<tol])
    loc_Hout = 0.0
    if any(xa > L-tol):
        loc_Hout = max(za[xa > L-tol])
    bs = mesh.comm.allreduce(loc_bs, op=MPI.MIN)
    Hout = mesh.comm.allreduce(loc_Hout, op=MPI.MAX)
    isslab = abs(bs) < tol and abs(Hin - Hout) < tol
    return (bs,Hout,isslab)

def processopts():
    import argparse
    parser = argparse.ArgumentParser(description='Generate .geo geometry-description file, suitable for meshing by Gmsh, for the outline of a glacier flowing over a step.  Also generates slab-on-slope geometry with -bs 0.0.')
    parser.add_argument('-o', metavar='FILE.geo', default='glacier.geo',
                        help='output file name (ends in .geo; default=glacier.geo)')
    parser.add_argument('-bs', type=float, default=120.0, metavar='X',
                        help='height of bed step (default=100 m)')
    parser.add_argument('-hmesh', type=float, default=100.0, metavar='X',
                        help='default target mesh spacing (default=100 m)')
    parser.add_argument('-refine', type=float, default=1.0, metavar='X',
                        help='refine resolution by this factor (default=1)')
    parser.add_argument('-refine_corner', type=float, default=4.0, metavar='X',
                        help='further local refinement at interior corner by this factor (default=4)')
    return parser.parse_args()

if __name__ == "__main__":
    import sys
    from datetime import datetime
    import platform

    args = processopts()
    commandline = " ".join(sys.argv[:])  # for save in comment in generated .geo
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    print('writing bedrock-step geometry to file %s ...' % args.o)
    geo = open(args.o, 'w')
    # header which records creation info
    geo.write('// geometry-description file created %s by %s using command\n//   %s\n\n'
              % (now,platform.node(),commandline) )
    # set "characteristic lengths" which are used by gmsh to generate triangles
    lc = args.hmesh / args.refine
    print('setting target mesh size of %g m' % lc)
    geo.write('lc = %f;\n' % lc)
    if abs(args.bs) >= 1.0:
        lc_corner = lc / args.refine_corner
        print('setting target mesh size of %g m at interior corner' % lc_corner)
        geo.write('lc_corner = %f;\n' % lc_corner)
    else:
        print('WARNING: bed step height |bs| < 1.0 so bed step is removed')
    # the rest
    writegeometry(geo,args.bs)
    geo.close()
    print('done')

