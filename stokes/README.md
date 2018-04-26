stokes/
=======

Copyright 2018  Ed Bueler

The codes and documents in this directory use (or assume the use of) [Firedrake](https://www.firedrakeproject.org/), [Gmsh](http://gmsh.info/), [PETSc](http://www.mcs.anl.gov/petsc/), and [Paraview](https://www.paraview.org/).  This material is both more-advanced and more-experimental than codes in `mfiles/`

Installation
------------

Install Gmsh and Paraview as desired, for instance by installing Debian or OSX packages.  Then follow the instructions at the [Firedrake download page](https://www.firedrakeproject.org/download.html) to install it.  In a normal installation one should use the PETSc which is installed by Firedrake, so there is no separate PETSc installation process needed.

Default usage
-------------

        $ ./genstepmesh.py glacier.geo             # create domain geometry
        $ gmsh -2 glacier.geo                      # mesh domain
        $ source ~/firedrake/bin/activate          # start Firedrake
        (firedrake) $ ./flowstep.py glacier.msh    # solve Stokes problem
        (firedrake) $ paraview glacier.pvd         # visualize

Note that script `genstepmesh.py` allows uniform refinement (`-refine X`) and additional refinement at the interior corner created by the bedrock step (`-refine_corner X`).

Visualization in Paraview is made easier by loading the state file `flowstep.pvsm`.

Slab-on-slope usage
-------------------

Set the height of the bedrock step to zero when creating the domain geometry:

        $ ./genstepmesh.py -bs 0.0 slab.geo
        $ gmsh -2 slab.geo
        $ source ~/firedrake/bin/activate
        (firedrake) $ ./flowstep.py -bs 0.0 -f slab

Firedrake info
--------------

  * `unset PETSC_DIR` and `unset PETSC_ARCH` may be needed before running `activate` when starting Firedrake
  * do `python3 firedrake-status` in the `firedrake/bin/` directory to see the current configuration of your firedrake installation

