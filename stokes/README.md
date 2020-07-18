stokes/
=======

Copyright 2018--2020 Ed Bueler

The programs in this directory use [Firedrake](https://www.firedrakeproject.org/), [Gmsh](http://gmsh.info/), [PETSc](http://www.mcs.anl.gov/petsc/), and [Paraview](https://www.paraview.org/).  They are more advanced, and more experimental, than the Matlab/Octave programs in [mfiles/](../mfiles/).  The solver is documented in [doc/](doc/).

Installation
------------

  * Install [Gmsh](http://gmsh.info/) and [Paraview](https://www.paraview.org/),
    for instance by installing Debian or OSX packages.
  * Follow the instructions at the
    [Firedrake download page](https://www.firedrakeproject.org/download.html)
    to install Firedrake.
  * Most users will only need the PETSc which is installed by Firedrake; no
    separate PETSc installation is needed.

Stokes-only usage
-----------------

This is the basic mode of solving the momentum conservation problem, that is, the Glen-Stokes system for the velocity and pressure.  The glacier surface is fixed.

Start [Firedrake](https://www.firedrakeproject.org/):

        $ source ~/firedrake/bin/activate

Generate the domain geometry using [domain.py](domain.py) and then mesh it using [Gmsh](http://gmsh.info/):

        (firedrake) $ ./domain.py -o glacier.geo      # create domain outline
        (firedrake) $ gmsh -2 glacier.geo             # writes glacier.msh

Run the solver [flow.py](flow.py) on the mesh, which will write velocity and pressure fields into `glacier.pvd`.

        (firedrake) $ ./flow.py -mesh glacier.msh

Now visualize using [Paraview](https://www.paraview.org/):

        (firedrake) $ paraview glacier.pvd

An alternate visualization plots surface values into an image file:

        (firedrake) $ ./flow.py -mesh glacier.msh -osurface surf.png

Slab-on-slope usage
-------------------

Set the height of the bedrock step to zero when creating the domain geometry:

        (firedrake) $ ./domain.py -bs 0.0 -o slab.geo
        (firedrake) $ gmsh -2 slab.geo
        (firedrake) $ ./flow.py -mesh slab.msh

In this mode the numerical error is displayed because the exact solution is known.  (See [notes/](../notes/) and/or [doc/](doc/) for the slab-on-slope solution.)

Surface evolution usage
-----------------------

By setting `-deltat` to a positive value, and choosing the number of time steps by `-m`, the surface will evolve according to the surface kinematical equation, in the case of zero mass balance, using explicit time stepping.  The time-stepping is not adaptive and it is up to the user to find time steps which are stable.

For a stable example,

        (firedrake) $ ./flow.py -mesh glacier.msh -deltat 10.0 -m 60

This writes variables (velocity,pressure,vertical\_displacement) into `glacier.pvd` at each time step, and [Paraview](https://www.paraview.org/) can show an animation.  The alternate visualization will plot final-time-only surface values of (h,u,w,h_t) into the image file:

        (firedrake) $ ./flow.py -mesh glacier.msh -deltat 10.0 -m 60 -osurface foo.png

Unstable examples are in the script [study/genunstable.sh](study/genunstable.sh).

Mesh refinement
---------------

The default `glacier.msh` mesh above has a typical mesh size of 100 m with grid resolution a factor of 4 of two finer near the interior corners created by the bedrock step, giving 25 m resolution at these corners.  Remember you can visualized these mesh files with [Gmsh](http://gmsh.info/): `gmsh glacier.msh`.

There are four refinement methods to get finer resolution:

1. One may use the script `domain.py` to create a finer or locally-finer mesh, before it is read by `flow.py`.  The script allows setting the target mesh size (`-hmesh H`) and/or setting a uniform refinement factor (`-refine X`).  Another option controls the factor used for additional refinement at the interior corner (`-refine_corner Y`).  (_The default case corresponds to_ `-hmesh 100 -refine 1 -refine_corner 4`.)  For example the following creates a mesh with target mesh size varying from 25 m to about 3 m near the interior corner.  The resulting grid has about 15 times as many elements as the default mesh:

        (firedrake) $ ./domain.py -refine 4 -refine_corner 8 -o fine1.geo
        (firedrake) $ gmsh -2 fine1.geo
        (firedrake) $ ./flow.py -mesh fine1.msh

2. One may use [Gmsh](http://gmsh.info/) to refine an existing mesh in `.msh`, specifically by splitting each triangular cell (element) into four similar triangles:

        (firedrake) $ ./domain.py -refine 2 -refine_corner 4 -o start.geo
        (firedrake) $ gmsh -2 start.geo
        (firedrake) $ gmsh -refine start.msh -o fine2.msh
        (firedrake) $ ./flow.py -mesh fine2.msh

    As with method 1, this refinement is done before the mesh is read by `flow.py`.

3. One may read an initial mesh and ask the `flow.py` script, i.e. ask [Firedrake](https://www.firedrakeproject.org/), to refine:

        (firedrake) $ ./flow.py -mesh start.msh -refine 1

    This is slightly-simpler usage compared to method 2, and with the same result.

4. For Stokes-only computations one may use _grid-sequencing_.  This solves the problem on the initially-read mesh, interpolates onto the next finer one, and then solves there

        (firedrake) $ ./flow.py -mesh start.msh -sequence 1

    This is slightly faster than methods 2 and 3, but the result should be the same.  A convergence test using this refinement method and the exact slab-on-a-slope solution is in [study/convergeslab.sh](study/convergeslab.sh).

Coupled steady-state usage
--------------------------

Using any of methods 1, 2, and 3, one can generate high-quality (i.e. 1 mm/10 day) steady states of the coupled momentum (Stokes) and mass conservation (surface-kinematical) equations:

        (firedrake) $ ./flow.py -mesh fine1.msh -deltat 10.0 -m 100 -s_snes_converged_reason
        (firedrake) $ ./flow.py -mesh fine2.msh -deltat 10.0 -m 100 -s_snes_converged_reason
        (firedrake) $ ./flow.py -mesh start.msh -refine 1 -deltat 10.0 -m 100 -s_snes_converged_reason

Getting help
------------

There are several ways to get help, but note that the `-mesh` option may be
required to get going:

  * `./flow.py -flowhelp` shows options specific to `flow.py` itself,
  * `./flow.py -mesh X.msh -help` lists a large number of PETSc options, and
  * `./flow.py -mesh X.msh -help intro` shows the PETSc version number.

Solver performance information
------------------------------

[Firedrake](https://www.firedrakeproject.org/) calls [PETSc](http://www.mcs.anl.gov/petsc/) to solve the Glen-Stokes equations.  Because this is a nonlinear problem, the [SNES object](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/index.html) from [PETSc](http://www.mcs.anl.gov/petsc/) is used.  Also, the solver is referred-to using option prefix `-s_`.  Therefore basic information on solver performance comes, for example, from asking for the number of iterations in each solver, for example

        (firedrake) $ ./flow.py -mesh glacier.msh -deltat 10.0 -m 10 -s_snes_monitor -s_snes_converged_reason

For another example, a direct linear solver for each Newton step, with a print-out of [PETSc](http://www.mcs.anl.gov/petsc/) information about the solver, is chosen by

        (firedrake) $ ./flow.py -mesh glacier.msh -s_snes_view -s_mat_type aij -s_ksp_type preonly -s_pc_type lu

Note the usual matrix type is `nest`, but it does not work directly with the LU solver.  Inside the [SNES](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/index.html) solver are [KSP (Krylov space)](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/index.html) and [PC (preconditioner)](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/index.html) components.

The solvers are packaged so that users have guidance on PETSc options.  The default is `-package SchurDirect`.  For another example, the best solver so far for high-resolution grids uses Schur decomposition and geometric multigrid:

        (firedrake) $ ./flow.py -mesh glacier.msh -sequence 4 -package SchurGMGSelfp -s_snes_converged_reason -s_ksp_converged_reason

For easy cases (i.e. grids with a smallish number of nodes) one can write-out the system matrix as a Matlab file `matrix.m` as follows:

        (firedrake) $ ./flow.py -mesh glacier.msh -s_mat_type aij -s_ksp_view_mat :matrix.m:ascii_matlab

See the [PETSc](http://www.mcs.anl.gov/petsc/) solver options in [momentummodel.py](momentummodel.py).

In time-stepping mode there is a second solver which computes the mesh vertical displacement field by solving Laplace's equation.  It has option prefix `-vd_` and it is defined, with default options, in [meshmotion.py](meshmotion.py).

For more on [PETSc](http://www.mcs.anl.gov/petsc/) and [Firedrake](https://www.firedrakeproject.org/) solvers see their online documentation, or see my book [_PETSc for PDEs_](https://github.com/bueler/p4pdes).

Info on installing Firedrake
-----------------------

  * You may need to `unset PETSC_DIR` and `unset PETSC_ARCH` before running `activate` when starting [Firedrake](https://www.firedrakeproject.org/).
  * Do `python3 firedrake-status` in the `firedrake/bin/` directory, after running `activate`, to see the current configuration of your [Firedrake](https://www.firedrakeproject.org/) installation.

