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

import numpy
import argparse
import sys
from datetime import datetime
import platform

commandline = " ".join(sys.argv[:])  # for save in comment in generated .geo
now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

# process options
parser = argparse.ArgumentParser(description='Generate .geo geometry-description file, suitable for meshing by Gmsh, for the outline of a glacier flowing over a step.  Also generates slab-on-slope geometry (use -bs 0.0).')
parser.add_argument('filename',
                    help='output file name (ends in .geo)')
parser.add_argument('-bs', type=float, default=120.0, metavar='X',
                    help='height of bed step (m)')
parser.add_argument('-hmesh', type=float, default=100.0, metavar='X',
                    help='default target mesh spacing')
parser.add_argument('-refine', type=float, default=1.0, metavar='X',
                    help='refine resolution by this factor')
parser.add_argument('-refine_corner', type=float, default=4.0, metavar='X',
                    help='locally refine at interior corner by this factor')
args = parser.parse_args()
bs = args.bs

# open .geo file and put in header which records creation info
geo = open(args.filename, 'w')
print('writing bedrock-step geometry to file %s ...' % args.filename)
geo.write('// geometry-description file created %s by %s using command\n//   %s\n\n'
          % (now,platform.node(),commandline) )

# domain distances in meters
H = 400.0        # input and output thickness (z)
L = 2000.0       # total along-flow length (x)
Ls = 1300.0      # location of bedrock step (x)
Ks = 150.0       # distance over which transition in surface happens (x)

# set "characteristic lengths" which are used by gmsh to generate triangles
lc = args.hmesh / args.refine
print('setting target mesh size of %g m' % lc)
geo.write('lc = %f;\n' % lc)
if abs(bs) >= 1.0:
    lc_corner = lc / args.refine_corner
    print('setting target mesh size of %g m at interior corner' % lc_corner)
    geo.write('lc_corner = %f;\n' % lc_corner)
else:
    print('WARNING: bed step height |bs| < 1.0 so bed step is removed')

# points on boundary
geo.write('Point(1) = {%f,%f,0,lc};\n' % (L,0.0))
geo.write('Point(2) = {%f,%f,0,lc};\n' % (L,H))
geo.write('Point(3) = {%f,%f,0,lc};\n' % (Ls+Ks,H))
# FIXME  between Point(3) and Point(4) could be a nice smooth curve
geo.write('Point(4) = {%f,%f,0,lc};\n' % (Ls-4*Ks,H+bs))
geo.write('Point(5) = {%f,%f,0,lc};\n' % (0.0,H+bs))
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
geo.write('Physical Line(41) = {11};\n')        # outflow (right side)
geo.write('Physical Line(42) = {12,13,14};\n')  # top (stress-free)
geo.write('Physical Line(43) = {15};\n')        # inflow (left side)
if abs(bs) >= 1.0:                              # base (no slip)
    geo.write('Physical Line(44) = {16,17,18};\n')
else:
    geo.write('Physical Line(44) = {16};\n')
geo.write('Physical Surface(51) = {31};\n')   # ensure all interior elements written(?)

# finish up
geo.close()
print('done')

