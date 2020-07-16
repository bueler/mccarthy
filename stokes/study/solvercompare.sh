#!/bin/bash
set -e

# compare solver performance on sequence of refining meshes

# notes:
#   1. Direct does not scale for higher res?
#   2. Testing of SchurGMGSelfp up to sequence 5 suggests optimality;
#      only slower than SchurDirect up to -sequence 3.

SEQ=2

../gendomain.py -o glacier.geo | grep writing
gmsh -2 glacier.geo | grep Writing
echo

for PAC in SchurDirect SchurGMGSelfp Direct; do
    CMD="../flow.py -mesh glacier.msh -sequence $SEQ -package $PAC -s_snes_converged_reason"
    echo $CMD

    rm -f tmp.txt
    /usr/bin/time -f "real %e" $CMD &> tmp.txt
    grep "grid-sequencing" tmp.txt
    grep "solve converged" tmp.txt
    grep "real" tmp.txt
done

