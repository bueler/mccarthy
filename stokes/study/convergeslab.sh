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
  CMD="mpiexec -n $P ../flow.py -mesh $1 -eps $2 -s_snes_converged_reason -s_snes_max_it 200"
  echo $CMD
  rm -f tmp.txt
  #/usr/bin/time -f "real %e" $CMD &> tmp.txt
  $CMD &> tmp.txt
  grep 'converged due to' tmp.txt
  grep 'flow speed' tmp.txt
  grep 'numerical errors' tmp.txt
  #grep real tmp.txt
  rm -f tmp.txt
}

EPS=0.0001
for REFINE in 0 1 2 4; do
    MESH=slab$REFINE
    if [ "$REFINE" -eq "0" ]; then
        ../gendomain.py -hmesh 160.0 -bs 0.0 -o $MESH.geo | grep setting
    else
        ../gendomain.py -refine $REFINE -bs 0.0 -o $MESH.geo | grep setting
    fi
    gmsh -2 $MESH.geo | grep nodes
    runcase $MESH.msh $EPS 
done

