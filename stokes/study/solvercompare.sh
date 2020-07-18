#!/bin/bash
set -e

# compare solver performance on sequence of refining meshes

# notes:
#   1. Direct does not converge for higher res
#   2. Testing of SchurGMGSelfp up to sequence 5 suggests optimality;
#      only slower than SchurDirect up to -sequence 3.

SEQ=3
SOLVERS="SchurDirect SchurGMGSelfp SchurGMGMass SchurGMGMassFull"

../gendomain.py -o glacier.geo | grep writing
gmsh -2 glacier.geo | grep Writing
echo

for PAC in $SOLVERS; do
    CMD="../flow.py -mesh glacier.msh -sequence $SEQ -package $PAC -s_snes_converged_reason -s_ksp_converged_reason"
    echo $CMD

    rm -f tmp.txt
    /usr/bin/time -f "real %e" $CMD &> tmp.txt
    grep "grid-sequencing" tmp.txt
    grep -B 1 "Nonlinear s_ solve converged" tmp.txt
    grep "real" tmp.txt
done

