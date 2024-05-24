mccarthy
========

Copyright 2010--2024  Ed Bueler

This repository contains slides, notes, and computer programs about numerical glacier and ice sheet modeling.  They have been used for the [International Summer School in Glaciology](http://glaciers.gi.alaska.edu/courses/summerschool), McCarthy, AK in years 2010, 2012, 2014, 2016, 2018, and 2022.


UNDER CONSTRUCTION!
-------------------

This content will be completely revised for the June 2024 school.  The current state is messy ...


slides and notes
----------------

The PDF slides [slides/slides-2024.pdf](slides/slides-2024.pdf) are new this year.  I believe they cover the essential material in a better way.  They are supported by Python codes in [py/](py/).

The older material is in PDF notes [notes/notes-2024.pdf](notes/notes-2024.pdf), previous slides [slides/slides-2022.pdf](slides/slides-2022.pdf), and Matlab/Octave programs [mfiles/](mfiles/).  

Python programs
---------------

The codes in subdirectory [py/](py/) solve surface kinematic equation and Stokes problems.  Download them either by cloning this repository or by getting a `.zip` or `.tar.gz` archive at the [releases page](https://github.com/bueler/mccarthy/releases).  For now the [new slides](slides/slides-2024.pdf) are the documentation, along with line comments.

Matlab/Octave programs (deprecated)
-----------------------------------

The codes in subdirectory [mfiles/](mfiles/) solve SIA and SSA problems.  They should work in either Matlab or [Octave](https://www.gnu.org/software/octave/); if not please report a bug, either by email or by using the [issues](https://github.com/bueler/mccarthy/issues) for this repository.  Download them either by cloning this repository or by getting a `.zip` or `.tar.gz` archive at the [releases page](https://github.com/bueler/mccarthy/releases).  The older notes and slides document these programs, but the programs also have help files (i.e. leading comments).  You are encouraged to actually run and modify them!

stokes solver
-------------

The Python tools in [stokes/](stokes/) are primarily for projects.  They solve a free-surface 2D Glen-Stokes flow over a bedrock step.  The workflow uses the following tools: [Firedrake](https://www.firedrakeproject.org/) (a finite element library), [Gmsh](http://gmsh.info/) (a mesh generator), [PETSc](http://www.mcs.anl.gov/petsc/) (a solver library), and [Paraview](https://www.paraview.org/) (for visualization).  See [stokes/README.md](stokes/README.md) for more information.


ancient versions
----------------

Much older versions (2009, 2010, 2012, 2014) lived in the repo https://github.com/bueler/karthaus
