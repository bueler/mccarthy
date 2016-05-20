#!/usr/bin/env python
# (C) 2015 Ed Bueler

import argparse
import sys
import numpy as np

import PetscBinaryIO as pbio
import petsc_conf

parser = argparse.ArgumentParser(description='Generate PETSc binary files with Gaussian-process generate random bedrock topograpy.')
parser.add_argument('outname', metavar='OUTNAME',
                    help='output PETSc binary file with variables x,y,topg')
parser.add_argument("--fixme", action="store_true",
                    help="fixme")
args = parser.parse_args()


FIXME FROM HERE: implement bedrand2d in python


# convert to PETSc-type vecs
xvec = x.view(pbio.Vec)
yvec = y.view(pbio.Vec)
topgvec = topg.flatten().view(pbio.Vec)

# open petsc binary file
io = pbio.PetscBinaryIO()

# write fields **in a particular order**  (names do not matter)
print "writing vars x,y,topg into %s ..." % args.outname
io.writeBinaryFile(args.outname, [xvec,yvec,topgvec,])

