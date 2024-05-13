#!/usr/bin/env python3
# (C) 2018--2024 Ed Bueler

# The purpose of this code is to generate a flow-line domain,
# suitable for meshing, on which we can solve the Stokes equations.
# See README.md and flow.py for usage in the Stokes context.
# Basic usage:
#   $ ./domain.py glacier.geo
#   $ gmsh -2 glacier.geo
# This generates glacier.msh.  You can use gmsh to inspect the mesh,
# or you can load it into python using firedrake.  For example:
#   $ source firedrake/bin/activate
#   (firedrake) $ ipython3
#   In [1]: from firedrake import *
#   In [2]: Mesh('glacier.msh')

# numbering of parts of boundary
bdryids = {'outflow' : 41,
           'top'     : 42,
           'inflow'  : 43,
           'base'    : 44}

def writegeometry(geo,bs, L, Lup, Ldown, Hin, Hout):
    # points on boundary, with target mesh densities
    Lmid = 0.5 * (Lup + Ldown)
    geo.write('Point(1) = {%f,%f,0,lc};\n' % (L,0.0))
    geo.write('Point(2) = {%f,%f,0,lc};\n' % (L,Hout))
    geo.write('Point(3) = {%f,%f,0,lc};\n' % (0.0,Hin))
    geo.write('Point(4) = {%f,%f,0,lc};\n' % (0.0,0.0))
    # if there is a bedrock step, *refine* in the interior corners
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

    # ensure all interior elements written
    geo.write('Physical Surface(51) = {31};\n')

def processopts():
    import argparse
    parser = argparse.ArgumentParser(description=
    '''Generate .geo geometry-description file, suitable for meshing by Gmsh,
    for the outline of a glacier flow domain with bedrock steps.  Also
    generates slab-on-slope geometry with -bs 0.0.
    ''')
    parser.add_argument('-bs', type=float, default=100.0, metavar='X',
                        help='height of bed step (default=100 m)')
    parser.add_argument('-Hin', type=float, default=400.0, metavar='X',
                        help='upstream thickness of ice (default=400 m)')
    parser.add_argument('-Hout', type=float, default=400.0, metavar='X',
                        help='downstream thickness of ice (default=400 m)')
    parser.add_argument('-hmesh', type=float, default=80.0, metavar='X',
                        help='default target mesh spacing (default=80 m)')
    parser.add_argument('-L', type=float, default=3000.0, metavar='X',
                        help='flow line length (default=3000 m)')
    parser.add_argument('-Lup', type=float, default=1500.0, metavar='X',
                        help='start of bedrock step (default=1500 m)')
    parser.add_argument('-Ldown', type=float, default=2000.0, metavar='X',
                        help='end of bedrock step (default=2000 m)')
    parser.add_argument('-o', metavar='FILE.geo', default='glacier.geo',
                        help='output file name (ends in .geo; default=glacier.geo)')
    parser.add_argument('-refine', type=float, default=1.0, metavar='X',
                        help='refine resolution by this factor (default=1)')
    parser.add_argument('-refine_corner', type=float, default=4.0, metavar='X',
                        help='further local refinement at interior corner by this factor (default=4)')
    return parser.parse_args()

if __name__ == "__main__":
    from datetime import datetime
    import sys, platform, subprocess

    args = processopts()
    commandline = " ".join(sys.argv[:])  # for save in comment in generated .geo
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    print('writing domain geometry to file %s ...' % args.o)
    geo = open(args.o, 'w')
    # header which records creation info
    geo.write('// geometry-description file created %s by %s\n'
              % (now,platform.node()) )
    geo.write('// command used:\n//   %s\n\n' % commandline)
    # set "characteristic lengths" which are used by gmsh to generate triangles
    lc = args.hmesh / args.refine
    geo.write('lc = %f;\n' % lc)
    if abs(args.bs) > 1.0:
        lc_corner = lc / args.refine_corner
    else:
        lc_corner = lc
    geo.write('lc_corner = %f;\n' % lc_corner)
    # write the rest of the .geo file
    writegeometry(geo, args.bs, args.L, args.Lup, args.Ldown, args.Hin, args.Hout)
    geo.close()
