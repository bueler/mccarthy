#!/bin/bash
set -e

# This script deliberately shows unstable evolution of the surface using
# explicit updating.  The first run uses default coarsening in the upper left
# corner (-coarsen_upperleft 2.0), and the resulting instability is initiated
# over the upstream end of the bed step.  The second run *refines* in the upper
# left corner and the instability is initated in that upper left corner.

# View the results with paraview and animation:
#   (firedrake) $ ./genunstable.sh
#   ...
#   (firedrake) $ paraview unstableA.pvd
#   (firedrake) $ paraview unstableB.pvd

../gendomain.py -o unstableA.geo
gmsh -2 unstableA.geo
../flow.py -mesh unstableA.msh -deltat 60.0 -m 20 -s_snes_converged_reason

../gendomain.py -coarsen_upperleft 0.5 -o unstableB.geo
gmsh -2 unstableB.geo
../flow.py -mesh unstableB.msh -deltat 40.0 -m 30 -s_snes_converged_reason

