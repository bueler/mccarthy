#!/usr/bin/env python3
# (C) 2022 Ed Bueler

# For usage see README.md.

# computational domain dimensions in meters
H0 = 400.0       # input (and initial output) thickness (z)
L = 3000.0       # total along-flow length (x)

# numbering of parts of boundary
bdryids = {'outflow' : 41,  # not used in this geometry
           'top'     : 42,
           'inflow'  : 43,
           'base'    : 44}

# define margin top shape as a parameterized curve for 0 <= t <= 1
def curve(t):
    return L * (1.0 - t**1.5), H0 * t

def writegeometry(geo,ds):
    assert (ds < L / 3.0)
    import numpy as np
    # points on boundary, with target mesh densities lc
    geo.write('Point(1) = {%f,%f,0,lc};\n' % (0.0,H0))
    geo.write('Point(2) = {%f,%f,0,lc};\n' % (0.0,0.0))
    geo.write('Point(3) = {%f,%f,0,lc};\n' % (L,0.0))
    # points on top defined by arclength scheme
    t, n = 0, 3
    x, y = curve(t)
    C, D = 8.0 * H0**3 / (27.0 * L**2), 9.0 * L**2 / (4.0 * H0**2)
    while x > L / 2.0:
        # from arclength spacing ds compute parameter spacing h
        tmp = ds / C + (1.0 + D * t)**1.5
        h = (tmp**(2.0/3.0) - 1.0) / D - t
        # update t, x, y and put point in file
        t = t + h
        x, y = curve(t)
        n = n + 1
        geo.write('Point(%d) = {%f,%f,0,lc};\n' % (n,x,y))
    # remaining points on top, back to y-axis
    m = int(np.ceil(x / ds))
    deltax = x / m
    for j in range(m-1):
        x = x - deltax
        y = H0 * (1.0 - x / L)**(2/3)
        n = n + 1
        geo.write('Point(%d) = {%f,%f,0,lc};\n' % (n,x,y))
    # lines along boundary
    for k in range(n-1):
       geo.write('Line(%d) = {%d,%d};\n' % (k+1,k+1,k+2))
    geo.write('Line(%d) = {%d,%d};\n' % (n,n,1))
    # line loop and surface allows defining a 2D mesh
    LineLoopN = n+1
    PlaneSurfaceN = n+2
    geo.write('Line Loop(%d) = {1' % LineLoopN)
    for k in range(n-1):
       geo.write(',%d' % (k+2))
    geo.write('};\n')
    geo.write('Plane Surface(%d) = {%d};\n' % (PlaneSurfaceN,LineLoopN))
    # "Physical" for marking boundary conditions
    geo.write('Physical Line(%d) = {1};\n' % bdryids['inflow'])
    geo.write('Physical Line(%d) = {2};\n' % bdryids['base'])
    geo.write('Physical Line(%d) = {3' % bdryids['top'])
    for k in range(4,n+1):
       geo.write(',%d' % k)
    geo.write('};\n')
    # ensure all interior elements written ... NEEDED!
    geo.write('Physical Surface(51) = {%d};\n' % PlaneSurfaceN)

def processopts():
    import argparse
    parser = argparse.ArgumentParser(description=
    '''Generate .geo geometry-description file, suitable for meshing by Gmsh.
    ''')
    parser.add_argument('-o', metavar='FILE.geo', default='margin.geo',
                        help='output file name (ends in .geo; default=margin.geo)')
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

    print('writing profile geometry to file %s ...' % args.o)
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
    writegeometry(geo,lc)
    geo.close()
    if args.testspew:
        result = subprocess.run(['cat', args.o], stdout=subprocess.PIPE)
        print(result.stdout)
