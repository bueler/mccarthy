#!/usr/bin/env python
# (C) 2015 Ed Bueler

import argparse
import sys

import numpy as np
from scipy.sparse.linalg import eigs
from numpy.random import randn

import PetscBinaryIO as pbio
import petsc_conf

parser = argparse.ArgumentParser(description=
'''Generate PETSc binary files with random (Gaussian-process-generated) periodic
bedrock topography.  Optionally show the result in a graphic.

Domain is rectangle square [0,Lx] x [0,Ly], with Nx x Ny grid.

Covariance function is
     K((x,y),(x',y')) = exp(- ( X^2 + Y^2 ) / (2*l^2))
where X = Lx * sin(pi(x-x')/Lx), Y = Ly * sin(pi(x-x')/Lx), and l is the length
scale.

The actual bedrock is shifted and scaled: if Z(x,y) is the sample of the
Gaussian process generated using the above covariance function then
     b(x,y) = A Z(x,y) + b0.
''',
#formatter_class=argparse.RawDescriptionHelpFormatter,  # irritation: want BOTH of these formatters
formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('outname', metavar='OUTNAME',
                    help='output PETSc binary file with variables x,y,b')
parser.add_argument("--amplitude", default=1000.0, metavar='A', type=float,
                    help="amplitude A of variation of b(x,y), in meters")
parser.add_argument("--lengthscale", default=20000.0, metavar='L', type=float,
                    help="distance scale l, in meters, used in periodized squared-exponential covariance function")
parser.add_argument("-Lx", default=100000.0, metavar='LX', type=float,
                    help="x extent of region, in meters")
parser.add_argument("-Ly", default=100000.0, metavar='LY', type=float,
                    help="y extent of region, in meters")
parser.add_argument("--meanelevation", default=2000.0, metavar='B0', type=float,
                    help="mean of b(x,y), in meters")
parser.add_argument("-Nx", default=20, metavar='NX', type=int,
                    help="number of grid points (and intervals) in x direction")
parser.add_argument("-Ny", default=20, metavar='NY', type=int,
                    help="number of grid points (and intervals) in y direction")
#FIXME: add ability to do many samples
#parser.add_argument("-P", default=1, metavar='P', type=int,
#                    help="number of beds to generate (samples)")
parser.add_argument("--plotit", action="store_true",
                    help="show with matplotlib and a wire frame")
args = parser.parse_args()

dx = args.Lx / float(args.Nx)
dy = args.Ly / float(args.Ny)
x = np.linspace(0,args.Lx-dx,args.Nx)
y = np.linspace(0,args.Ly-dy,args.Ny)

xx, yy = np.meshgrid(x,y)
xxx = xx.flatten()
yyy = yy.flatten()

N = args.Nx * args.Ny  # total number of points
Sigma = np.zeros((N,N))

l = args.lengthscale
for k in range(args.Ny):
    for j in range(args.Nx):
        SX = args.Lx * np.sin(np.pi * abs(xxx-x[j]) / args.Lx)
        SY = args.Ly * np.sin(np.pi * abs(yyy-y[k]) / args.Ly)
        Sigma[:][(k-1)*args.Nx + j] = np.exp(- (SX**2 + SY**2) / (2.0 * l * l) )

D, V = eigs(Sigma, k=200)
V = np.real(V)
Dhalf = np.sqrt(abs(D))
print 'Dhalf range is from %.6f to %.6f' % (min(Dhalf),max(Dhalf))

A = args.amplitude
b0 = args.meanelevation
Z = np.dot(V, Dhalf * np.dot(V.transpose(),randn(N)))
b = A * np.reshape(Z,(args.Nx,args.Ny)) + b0

if args.plotit:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from mpl_toolkits.mplot3d import axes3d
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(xx/1000.0,yy/1000.0,b,rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0)
    #plt.pcolor(xx/1000.0,yy/1000.0,b)
    ax.set_xlabel('x in km')
    ax.set_ylabel('y in km')
    ax.set_zlabel('b(x,y) in m')
    ax.set_zlim3d(0.0, max(b.flatten()))
    fig.colorbar(surf,shrink=0.8,aspect=8)
    plt.show()

# convert to PETSc-type vecs
xvec = x.view(pbio.Vec)
yvec = y.view(pbio.Vec)
bvec = b.flatten().view(pbio.Vec)

# open petsc binary file
io = pbio.PetscBinaryIO()

# write fields **in a particular order**  (names do not matter)
print "writing vars x,y,b into %s ..." % args.outname
io.writeBinaryFile(args.outname, [xvec,yvec,bvec,])

