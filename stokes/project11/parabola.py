#!/usr/bin/env python3
# (C) 2022 Ed Bueler

# Usage:
#   $ ./parabola.py margin.geo
#   $ gmsh -2 margin.geo
# which generates margin.msh.

# computational domain dimensions in meters
H0 = 400.0       # input (and initial output) thickness (z)
L = 3000.0       # total along-flow length (x)

# numbering of parts of boundary; these values match flow.py etc.
bdryids = {'outflow' : 41,  # not used in this geometry
           'top'     : 42,
           'inflow'  : 43,
           'base'    : 44}

def writegeometry(geo):
    # points on boundary, with target mesh densities
    geo.write('Point(1) = {%f,%f,0,lc};\n' % (0.0,0.0))
    geo.write('Point(2) = {%f,%f,0,lc};\n' % (L,0.0))
    FIXME geo.write('Point(2) = {%f,%f,0,lc};\n' % (L,Hout))

    # lines along boundary
    FIXME
    geo.write('Line(11) = {1,2};\n')
    geo.write('Line(12) = {2,3};\n')
    geo.write('Line(13) = {3,4};\n')
    geo.write('Line(14) = {4,1};\n')
    geo.write('Line Loop(21) = {11,12,13,14};\n')

    # surface allows defining a 2D mesh
    geo.write('Plane Surface(31) = {21};\n')

    # "Physical" for marking boundary conditions
    FIXME
    geo.write('Physical Line(%d) = {12};\n' % bdryids['top'])
    geo.write('Physical Line(%d) = {13};\n' % bdryids['inflow'])
    geo.write('Physical Line(%d) = {11};\n' % bdryids['base'])

    # ensure all interior elements written ... NEEDED!
    geo.write('Physical Surface(51) = {31};\n')

def processopts():
    import argparse
    parser = argparse.ArgumentParser(description=
    '''Generate .geo geometry-description file, suitable for meshing by Gmsh.
    ''')
    parser.add_argument('-o', metavar='FILE.geo', default='glacier.geo',
                        help='output file name (ends in .geo; default=glacier.geo)')
    parser.add_argument('-hmesh', type=float, default=80.0, metavar='X',
                        help='default target mesh spacing (default=80 m)')
    parser.add_argument('-testspew', action='store_true',
                        help='write .geo contents, w/o header, to stdout', default=False)  # just for testing
    return parser.parse_args()

if __name__ == "__main__":
    from datetime import datetime
    import sys, platform, subprocess

    args = processopts()
    commandline = " ".join(sys.argv[:])  # for save in comment in generated .geo
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    print('writing parabola geometry to file %s ...' % args.o)
    geo = open(args.o, 'w')
    # header which records creation info
    if not args.testspew:
        geo.write('// geometry-description file created %s by %s using command\n//   %s\n\n'
                  % (now,platform.node(),commandline) )
    # set "characteristic lengths" which are used by gmsh to generate triangles
    lc = args.hmesh
    print('setting target mesh size of %g m' % lc)
    geo.write('lc = %f;\n' % lc)

    # the rest
    writegeometry(geo)
    geo.close()
    if args.testspew:
        result = subprocess.run(['cat', args.o], stdout=subprocess.PIPE)
        print(result.stdout)
