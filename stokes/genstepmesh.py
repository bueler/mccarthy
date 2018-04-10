#!/usr/bin/env python3
# (C) 2018 Ed Bueler

# Usage:
#   $ ./genstepmesh.py glacier.geo
#   $ gmsh -2 glacier.geo
#   ...
# which generates glacier.msh.  (Use "gmsh glacier.geo" for GUI version.)
# Then load mesh in python/firedrake:
#   $ source firedrake/bin/activate
#   (firedrake) $ ipython
#   ...
#   In [1]: from firedrake import *
#   In [2]: Mesh('glacier.msh')
# Main purpose is to solve the Stokes equations on this domain.  See
# linflowstep.py:
#   $ source firedrake/bin/activate
#   (firedrake) $ ./linflowstep.py

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
parser.add_argument('-refine', type=float, default=1.0, metavar='X',
                    help='refine default resolution by this factor')
parser.add_argument('-refine_corner', type=float, default=4.0, metavar='X',
                    help='locally refine at interior corner by this factor')
args = parser.parse_args()
bs = args.bs

# open .geo file and put in header which records creation info
geo = open(args.filename, 'w')
print('writing file %s ...' % args.filename)
geo.write('// geometry-description file created %s by %s using command\n//   %s\n\n'
          % (now,platform.node(),commandline) )

# glacier domain parameters; distance in meters
H = 400.0        # input and output thickness (z)
L = 2000.0       # total along-flow length (x)
Ls = 1300.0      # location of bedrock step (x)
Ks = 150.0       # distance over which transition in surface happens (x)

# mesh parameters
lc = 100.0 / args.refine
geo.write('lc = %f;\n' % lc)    # put target characteristic mesh size

# points on boundary
geo.write('Point(1) = {%f,%f,0,lc};\n' % (Ls,0.0))
geo.write('Point(2) = {%f,%f,0,lc};\n' % (L,0.0))
geo.write('Point(3) = {%f,%f,0,lc};\n' % (L,H))
geo.write('Point(4) = {%f,%f,0,lc};\n' % (Ls+Ks,H))
# FIXME  between Point(4) and Point(5) should be a nice smooth curve ... but what?
geo.write('Point(5) = {%f,%f,0,lc};\n' % (Ls-4*Ks,H+bs))
geo.write('Point(6) = {%f,%f,0,lc};\n' % (0.0,H+bs))
geo.write('Point(7) = {%f,%f,0,lc};\n' % (0.0,bs))
if abs(bs) >= 1.0:
    lc_corner = lc / args.refine_corner
    geo.write('lc_corner = %f;\n' % lc_corner)
    geo.write('Point(8) = {%f,%f,0,lc_corner};\n' % (Ls,bs))
else:
    print('WARNING: bed step height |bs| < 1.0 so bed step removed')

# lines along boundary
loop = '9'
for j in range(6):
    geo.write('Line(%d) = {%d,%d};\n' % (j+9,j+1,j+2))
    if j > 0:
        loop += ',%d' % (j+9)
if abs(bs) >= 1.0:
    geo.write('Line(15) = {7,8};\n')
    geo.write('Line(16) = {8,1};\n')
    loop += ',15,16'
else:
    geo.write('Line(15) = {7,1};\n')
    loop += ',15'

# loop and surface allows defining a 2D mesh
geo.write('Line Loop(17) = {%s};\n' % loop)
geo.write('Plane Surface(21) = {17};\n')

# "Physical" for marking boundary conditions
geo.write('Physical Line(31) = {11,12,13};\n')  # stress-free boundary (top)
geo.write('Physical Line(33) = {14};\n')        # inflow (left side)
geo.write('Physical Line(34) = {10};\n')        # outflow (right side)
if abs(bs) >= 1.0:                              # base of ice
    geo.write('Physical Line(32) = {9,15,16};\n')
else:
    geo.write('Physical Line(32) = {9,15};\n')

# not clear following is needed
geo.write('Physical Surface(41) = {21};\n')   # ensure all interior elements written(?)

# finish up
geo.close()

