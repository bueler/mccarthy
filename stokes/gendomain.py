#!/usr/bin/env python3
# (C) 2018 Ed Bueler

# Usage:
#   $ ./gendomain.py glacier.geo
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
# README.md and flow.py.

# immutable domain distances in meters
Hin = 400.0      # input (and initial output) thickness (z)
L = 3000.0       # total along-flow length (x)
Lup = 1500.0     # location of bedrock step up (x)
Ldown = 2000.0   # location of bedrock step down (x);  Ldown > Lup

# numbering of parts of boundary
bdryids = {'outflow' : 41,
           'top'     : 42,
           'inflow'  : 43,
           'base'    : 44}

def writegeometry(geo,bs):
    # points on boundary
    Hout = Hin
    Lmid = 0.5 * (Lup + Ldown)
    geo.write('Point(1) = {%f,%f,0,lc};\n' % (L,0.0))
    geo.write('Point(2) = {%f,%f,0,lc};\n' % (L,Hout))
    geo.write('Point(3) = {%f,%f,0,lc};\n' % (0.0,Hin))
    geo.write('Point(4) = {%f,%f,0,lc};\n' % (0.0,0.0))
    if abs(bs) > 1.0:
        geo.write('Point(5) = {%f,%f,0,lc};\n' % (Lup,0.0))
        geo.write('Point(6) = {%f,%f,0,lc_corner};\n' % (Lup,bs))
        geo.write('Point(7) = {%f,%f,0,lc};\n' % (Lmid,bs))
        geo.write('Point(8) = {%f,%f,0,lc_corner};\n' % (Ldown,bs))
        geo.write('Point(9) = {%f,%f,0,lc};\n' % (Ldown,0.0))

    # lines along boundary
    geo.write('Line(11) = {1,2};\n')
    geo.write('Line(12) = {2,3};\n')
    geo.write('Line(13) = {3,4};\n')
    if abs(bs) > 1.0:
        geo.write('Line(14) = {4,5};\n')
        geo.write('Line(15) = {5,6};\n')
        geo.write('Line(16) = {6,7};\n')
        geo.write('Line(17) = {7,8};\n')
        geo.write('Line(18) = {8,9};\n')
        geo.write('Line(19) = {9,1};\n')
        geo.write('Line Loop(21) = {11,12,13,14,15,16,17,18,19};\n')
    else:
        geo.write('Line(14) = {4,1};\n')
        geo.write('Line Loop(21) = {11,12,13,14};\n')

    # surface allows defining a 2D mesh
    geo.write('Plane Surface(31) = {21};\n')

    # "Physical" for marking boundary conditions
    geo.write('Physical Line(%d) = {11};\n' % bdryids['outflow'])
    geo.write('Physical Line(%d) = {12};\n' % bdryids['top'])
    geo.write('Physical Line(%d) = {13};\n' % bdryids['inflow'])
    if abs(bs) > 1.0:
        geo.write('Physical Line(%d) = {14,15,16,17,18,19};\n' % bdryids['base'])
    else:
        geo.write('Physical Line(%d) = {14};\n' % bdryids['base'])

    # ensure all interior elements written ... NEEDED!
    geo.write('Physical Surface(51) = {31};\n')

# dynamically extract geometry making these definitions (tolerance=1cm):
#   bs      = height of bedrock step     = (min z-coordinate over Lup < x < Ldown)
#   Hout    = ice thickness at output    = (max z-coordinate at x=L)
# in parallel no process owns the whole mesh so MPI_Allreduce() is needed
def getdomaindims(mesh,tol=0.01):
    from mpi4py import MPI
    # mesh needs to be a Mesh from Firedrake
    xa = mesh.coordinates.dat.data_ro[:,0]  # .data_ro acts like VecGetArrayRead
    za = mesh.coordinates.dat.data_ro[:,1]
    loc_bs = 9.99e99
    xinmid = (xa > Lup) * (xa < Ldown)
    if any(xinmid):
        loc_bs = min(za[xinmid])
    bs = mesh.comm.allreduce(loc_bs, op=MPI.MIN)
    loc_Hout = 0.0
    if any(xa > L-tol):
        loc_Hout = max(za[xa > L-tol])
    Hout = mesh.comm.allreduce(loc_Hout, op=MPI.MAX)
    return (bs,Hout)

def processopts():
    import argparse
    parser = argparse.ArgumentParser(description=
    '''Generate .geo geometry-description file, suitable for meshing by Gmsh, for
    the outline of a glacier flow domain with bedrock steps.  Also generates
    slab-on-slope geometry with -bs 0.0.
    ''')
    parser.add_argument('-o', metavar='FILE.geo', default='glacier.geo',
                        help='output file name (ends in .geo; default=glacier.geo)')
    parser.add_argument('-bs', type=float, default=100.0, metavar='X',
                        help='height of bed step (default=100 m)')
    parser.add_argument('-hmesh', type=float, default=80.0, metavar='X',
                        help='default target mesh spacing (default=80 m)')
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

    print('writing domain geometry to file %s ...' % args.o)
    geo = open(args.o, 'w')
    # header which records creation info
    geo.write('// geometry-description file created %s by %s using command\n//   %s\n\n'
              % (now,platform.node(),commandline) )
    # set "characteristic lengths" which are used by gmsh to generate triangles
    lc = args.hmesh / args.refine
    print('setting target mesh size of %g m' % lc)
    geo.write('lc = %f;\n' % lc)
    if abs(args.bs) > 1.0:
        lc_corner = lc / args.refine_corner
        print('setting target mesh size of %g m at interior corners' % lc_corner)
    else:
        lc_corner = lc
    geo.write('lc_corner = %f;\n' % lc_corner)
    # the rest
    writegeometry(geo,args.bs)
    geo.close()

