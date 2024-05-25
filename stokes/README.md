# stokes/

Copyright 2018--2024 Ed Bueler

## TODO

  * add pytest tests
  * compare, and avoid unnecessary duplication with, [stokes-ice-tutorial](https://github.com/bueler/stokes-ice-tutorial)

## Introduction

The Python programs in this directory use advanced, open-source libraries and tools:

  * [Firedrake](https://www.firedrakeproject.org/), a Python finite element library
     * [PETSc](https://petsc.org/release/), a solver library typically installed by Firedrake
  * [Gmsh](http://gmsh.info/), a mesher
  * [Paraview](https://www.paraview.org/), for visualization of mesh functions

They are documented by the current README and by `doc/doc.pdf` (also in this directory).

## Installation

  * Install [Gmsh](http://gmsh.info/) and [Paraview](https://www.paraview.org/),
    for instance by installing Debian or OSX packages.
  * Follow the instructions at the
    [Firedrake download page](https://www.firedrakeproject.org/download.html)
    to install Firedrake.

## Usage

Activate your [Firedrake](https://www.firedrakeproject.org/) virtual environment (venv) first:

        $ source ~/firedrake/bin/activate

Generate the domain geometry using `domain.py` and then mesh it using [Gmsh](http://gmsh.info/):

        (firedrake) $ ./domain.py -o glacier.geo      # create domain outline
        (firedrake) $ gmsh -2 glacier.geo             # writes glacier.msh

You can visualize the generated mesh with [Gmsh](http://gmsh.info/):

        (firedrake) $ gmsh glacier.msh

Now run the Stokes solver `flow.py` on the mesh, which will write velocity and pressure fields into `glacier.pvd`.

        (firedrake) $ ./flow.py -mesh glacier.msh -o glacier.pvd

Now visualize the result using [Paraview](https://www.paraview.org/):

        (firedrake) $ paraview glacier.pvd

An alternate visualization plots surface values into an image file `surf.png`:

        (firedrake) $ ./flow.py -mesh glacier.msh -osurface surf.png

### Slab-on-slope usage

Set the height of the bedrock step to zero when creating the domain geometry:

        (firedrake) $ ./domain.py -bs 0.0 -o slab.geo
        (firedrake) $ gmsh -2 slab.geo
        (firedrake) $ ./flow.py -slab -mesh slab.msh

Here the numerical error is displayed because the exact solution is known.  See `doc/doc.pdf` for the slab-on-slope solution.

### Mesh refinement

The default mesh above (`glacier.msh`) has a typical mesh size of 80 m with grid resolution a factor of 4 finer near the interior corners created by the bedrock step, giving 20 m resolution at these corners.

There are four refinement methods to get finer resolution:

1. One may use the script `domain.py` to create a refined (or locally-refined) mesh, before it is read by `flow.py`.  `domain.py` allows setting the target mesh size (`-hmesh H`) and/or setting a uniform refinement factor (`-refine X`).  Another option controls the additional refinement at the interior corner (`-refine_corner Y`, where `Y` is a factor).  The default case corresponds to `-hmesh 80 -refine 1 -refine_corner 4`.  The following creates a refined mesh with target mesh size varying from 20 m to about 3 m near the interior corner.  The resulting grid has about 15 times as many elements as the default mesh:

        (firedrake) $ ./domain.py -refine 4 -refine_corner 8 -o fine1.geo
        (firedrake) $ gmsh -2 fine1.geo
        (firedrake) $ ./flow.py -mesh fine1.msh

2. One may use [Gmsh](http://gmsh.info/) to refine an existing mesh in `.msh`, specifically by splitting each triangular cell (element) into four similar triangles:

        (firedrake) $ ./domain.py -refine 2 -refine_corner 4 -o start.geo
        (firedrake) $ gmsh -2 start.geo
        (firedrake) $ gmsh -refine start.msh -o fine2.msh
        (firedrake) $ ./flow.py -mesh fine2.msh

    As with method 1, this refinement is done before the mesh is read by `flow.py`.

3. One may read an initial mesh and ask the `flow.py` script, via [Firedrake](https://www.firedrakeproject.org/), to refine:

        (firedrake) $ ./flow.py -mesh start.msh -refine 1

    This is slightly-simpler usage compared to method 2, and with the same result.

4. One may use _grid-sequencing_.  This solves the problem on the initially-read mesh, interpolates onto the next finer one, and then solves there

        (firedrake) $ ./flow.py -mesh start.msh -sequence 1

    This is slightly faster than methods 2 and 3, but the result should be the same.

### Solver performance information

[Firedrake](https://www.firedrakeproject.org/) calls [PETSc](http://www.mcs.anl.gov/petsc/) to solve the Glen-Stokes equations.  Because this is a nonlinear problem, a PETSc [SNES solver object](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/index.html) does the Newton iterations.  This solver is addressed using option prefix `-s_`.

Basic information on solver performance comes from reporting on the Newton iterations in each Stokes solve:

        (firedrake) $ ./flow.py -mesh glacier.msh -s_snes_monitor -s_snes_converged_reason

A detailed description of the PETSc solver comes from `-s_snes_view`.  Inside the [SNES](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/index.html) solver are [KSP (Krylov space)](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/index.html) and [PC (preconditioner)](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/index.html) components.

For easy cases (i.e. grids with a smallish number of nodes) one can write-out the system matrix as a Matlab file `foo.m` as follows:

        (firedrake) $ ./flow.py -mesh glacier.msh -s_ksp_view_mat :foo.m:ascii_matlab

Note `foo.m` is now a command in Matlab/Octave, and one may use `spy` to see the structure.  Or the matrix can be read into Python.

For more on [Firedrake](https://www.firedrakeproject.org/) and [PETSc](http://www.mcs.anl.gov/petsc/) solvers see their respective online documentation.  Or see my book!

## Getting help

There are several ways to get help:

  * `./flow.py -flowhelp` shows options specific to `flow.py` itself
  * `./flow.py -mesh glacier.msh -help` lists a large number of PETSc options
  * `./flow.py -mesh glacier.msh -help intro` shows the PETSc version number

Note that in the latter two usages the `-mesh` option is required in order to advance the program state until a point where a PETSc solver exists.

## Testing

Do this in the current directory (`mccarthy/stokes/`):

        $ pytest .
