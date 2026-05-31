mccarthy/py/
============

The Python codes in this subdirectory solve surface kinematic equation (SKE) and shallow ice approximation (SIA) problems.  The primary Python programs, namely [surface.py](surface.py) and [shallowuw.py](shallowuw.py) use only the numpy and matplotlib libraries.

The [notes](../notes/notes-2026.pdf) are the documentation for these programs, along with line comments in the `.py` source codes.  You are encouraged to actually run and modify these them!  Feel free to email me about how they work.

The [surface2d.py](surface2d.py) program uses [Firedrake](https://www.firedrakeproject.org/) (a finite element library) which calls [PETSc](http://www.mcs.anl.gov/petsc/) for equation solvers.

The codes in [stokes/](stokes/) also use Firedrake to solve a 2D Glen-law Stokes flow over a bedrock step.  The workflow for the Stokes codes also use [Gmsh](http://gmsh.info/) (a mesh generator) and [Paraview](https://www.paraview.org/) (for visualization).  See [stokes/doc.pdf](stokes/doc.pdf) for more information.  See also my [tutorial for using Firedrake to solve the Stokes equations](https://github.com/bueler/stokes-ice-tutorial), which is in a separate repository.

Please report any bugs in the codes either by email or by using the [issues](https://github.com/bueler/mccarthy/issues) for this repository.
