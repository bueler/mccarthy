#!/bin/bash
set -e

P=$1

# measure convergence on sequence of refining slab geometries
# note viscosity regularization relevant
# after setting-up firedrake in parent directory, run as:
#    ./convergeslab.sh P &> convergeslabP.txt
# where P is number of processes

function runcase() {
  CMD="mpiexec -n $P ../flowstep.py -eps $2 $1 -s_snes_max_it 200"
  echo $CMD
  rm -f tmp.txt
  #/usr/bin/time -f "real %e" $CMD &> tmp.txt
  $CMD &> tmp.txt
  grep 'converged due to' tmp.txt
  grep 'maximum velocity' tmp.txt
  grep 'numerical error' tmp.txt
  #grep real tmp.txt
  rm -f tmp.txt
}

EPS=0.0001
for REFINE in 0 1 2 4 8; do
    MESH=slab$REFINE
    if [ "$REFINE" -eq "0" ]; then
        ../genstepmesh.py -hmesh 200.0 -bs 0.0 -o $MESH.geo | grep setting
    else
        ../genstepmesh.py -refine $REFINE -bs 0.0 -o $MESH.geo | grep setting
    fi
    gmsh -2 $MESH.geo | grep vertices
    runcase $MESH.msh $EPS 
done

