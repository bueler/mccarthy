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

Generate the domain geometry and the mesh:

        $ ./gendomain.py -o glacier.geo            # creates domain outline
        $ gmsh -2 glacier.geo                      # meshes domain and generates glacier.msh

Start Firedrake and run the solver:

        $ source ~/firedrake/bin/activate
        (firedrake) $ ./flow.py glacier.msh        # solves Stokes problem for velocity and pressure

This writes variables (velocity,pressure) into `glacier.pvd`.  Now visualize:

        (firedrake) $ paraview glacier.pvd

Alternative visualization is to plot surface values of (h,u,w) into an image file:

        (firedrake) $ ./flow.py -osurface foo.png glacier.msh

Slab-on-slope usage
-------------------

Set the height of the bedrock step to zero when creating the domain geometry:

        (firedrake) $ ./gendomain.py -bs 0.0 -o slab.geo
        (firedrake) $ gmsh -2 slab.geo
        (firedrake) $ ./flow.py slab.msh

Surface evolution usage
-----------------------

By setting `-deltat` to a positive value, and choosing the number of time steps by `-m`, the surface will evolve.  For example,

        (firedrake) $ ./flow.py -deltat 10.0 -m 60 glacier.msh

This writes variables (velocity,pressure,vertical\_displacement) into `glacier.pvd` at each time step.  Paraview can show an animation.  The alternate visualization will plot surface values of (h,u,w,h_t) into the image file:

        (firedrake) $ ./flow.py -osurface foo.png glacier.msh

Mesh refinement
---------------

The default mesh has a typical mesh size of 100 m with refinement by a factor of 4 near the interior corner created by the bedrock step.  (Giving 25 m resolution at the corner.)

Options to script `gendomain.py` allow setting the target mesh size and setting a uniform refinement factor: `-hmesh H -refine X`.  Another option controls the additional refinement at the interior corner: `-refine_corner Y`.  The default case corresponds to `-hmesh 100 -refine 1 -refine_corner 4`.

For example the following creates a mesh with target mesh size varying from 25 m to about 3 m near the interior corner.  The resulting grid has about 15 times as many elements as the default mesh:

        (firedrake) $ ./gendomain.py -refine 4 -refine_corner 8 -o finer.geo
        (firedrake) $ gmsh -2 finer.geo
        (firedrake) $ ./flow.py finer.msh

Solver performance information
------------------------------

There are two PETSc solvers at each time step in surface evolution mode.  One, with prefix `s_`, solves the Stokes equations from the current geometry.  It is nonlinear and computes the velocity and pressure.  The other solver, with prefix `t_`, computes the vertical displacement of the mesh to solve the surface kinematical equation.  Basic information on solver performance comes from asking for the number of iterations in each solver, for example

        (firedrake) $ ./flow.py -deltat 10.0 -m 60 glacier.msh -s_snes_converged_reason -t_ksp_converged_reason

Note that all PETSc options must come _after_ the options to the script `flow.py`.  Also, the `t_` solver is not used in the default diagnostic-only mode (i.e. when there is no positive value supplied as `-deltat`).

Info on using Firedrake
-----------------------

  * `unset PETSC_DIR` and `unset PETSC_ARCH` may be needed before running `activate` when starting Firedrake
  * do `python3 firedrake-status` in the `firedrake/bin/` directory to see the current configuration of your firedrake installation

