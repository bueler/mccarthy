stokes/
=======

Copyright 2018 Ed Bueler

The example in this directory uses [Firedrake](https://www.firedrakeproject.org/), [Gmsh](http://gmsh.info/), [PETSc](http://www.mcs.anl.gov/petsc/), and [Paraview](https://www.paraview.org/).  It is more advanced, and more experimental, than codes in `mfiles/`

Installation
------------

  * Install [Gmsh](http://gmsh.info/) and [Paraview](https://www.paraview.org/),
    for instance by installing Debian or OSX packages.
  * Follow the instructions at the
    [Firedrake download page](https://www.firedrakeproject.org/download.html)
    to install Firedrake.
  * Most users will only need the PETSc which is installed by Firedrake; no
    separate PETSc installation is needed.

Default Stokes-only usage
-------------------------

        $ ./genstepmesh.py -o glacier.geo          # create domain geometry
        $ gmsh -2 glacier.geo                      # mesh domain
        $ source ~/firedrake/bin/activate          # start Firedrake
        (firedrake) $ ./flowstep.py glacier.msh    # solve Stokes problem
        (firedrake) $ paraview glacier.pvd         # visualize

Visualization in Paraview is made easier by loading the state file `flowstep.pvsm`.

Slab-on-slope usage
-------------------

Set the height of the bedrock step to zero when creating the domain geometry:

        (firedrake) $ ./genstepmesh.py -bs 0.0 -o slab.geo
        (firedrake) $ gmsh -2 slab.geo
        (firedrake) $ ./flowstep.py slab.msh

Surface evolution usage
-----------------------

FIXME

Mesh refinement
---------------

The default mesh has a typical mesh size of 100 m with refinement by a factor of 4 near the interior corner created by the bedrock step.  (Giving 25 m resolution at the corner.)

Options to script `genstepmesh.py` allow setting the mesh size and a uniform refinement factor: `-hmesh H -refine X`.  Another option controls the additional refinement at the interior corner: `-refine_corner Y`.  The default case corresponds to `-hmesh 100 -refine 1 -refine_corner 4`.

For example the following creates a mesh with target mesh size varying from 25 m to about 3 m near the interior corner.  The resulting grid has about 15 times as many elements as the default mesh:

        (firedrake) $ ./genstepmesh.py -refine 4 -refine_corner 8 -o finer.geo
        (firedrake) $ gmsh -2 finer.geo
        (firedrake) $ ./flowstep.py finer.msh

Firedrake info
--------------

  * `unset PETSC_DIR` and `unset PETSC_ARCH` may be needed before running `activate` when starting Firedrake
  * do `python3 firedrake-status` in the `firedrake/bin/` directory to see the current configuration of your firedrake installation

