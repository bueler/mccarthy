#!/bin/bash
set -e

# measure convergence relative to slab-on-slope on sequence of refining meshes

# after setting-up firedrake in parent directory, run as:
#    ./convergeslab.sh P &> convergeslabP.txt
# where P is number of processes

# FIXME multigrid should help here with reducing KSP iterations

# FIXME parallel runs are not really worthwhile with LU solvers ... redo after
#       solver testing

P=$1
EPS=0.00001  # viscosity regularization important to get close to exact
SEQUENCE=4


../gendomain.py -hmesh 160.0 -bs 0.0 -coarsen_upperleft 1.0 -o base.geo
gmsh -2 base.geo | grep Running

CMD="mpiexec -n $P ../flow.py -mesh base.msh -eps $EPS -sequence $SEQUENCE -s_snes_converged_reason -s_snes_max_it 200"
echo $CMD
$CMD

