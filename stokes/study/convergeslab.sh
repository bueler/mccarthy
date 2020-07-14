#!/bin/bash
set -e

P=$1

# measure convergence on sequence of refining slab geometries
# note viscosity regularization relevant: -eps ...
# after setting-up firedrake in parent directory, run as:
#    ./convergeslab.sh P &> convergeslabP.txt
# where P is number of processes

# grid-sequencing (and multigrid) would help here with reducing SNES
# iterations, but flow.py does not set up a grid hierarchy

function runcase() {
  CMD="mpiexec -n $P ../flow.py -mesh $1 -eps $2 -refine $3 -s_snes_converged_reason -s_snes_max_it 200"
  echo $CMD
  rm -f tmp.txt
  #/usr/bin/time -f "real %e" $CMD &> tmp.txt
  $CMD &> tmp.txt
  grep 'elements' tmp.txt
  grep 'converged due to' tmp.txt
  grep 'flow speed' tmp.txt
  grep 'numerical errors' tmp.txt
  #grep real tmp.txt
  rm -f tmp.txt
}

EPS=0.0001
../gendomain.py -hmesh 160.0 -bs 0.0 -coarsen_upperleft 1.0 -o base.geo
gmsh -2 base.geo | grep Running
for REFINE in 0 1 2 3; do
    runcase base.msh $EPS $REFINE
done

