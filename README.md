mccarthy
========

Copyright 2010--2023  Ed Bueler

This repository contains slides, notes, and computer programs about numerical glacier and ice sheet modeling.  They have been used for the [International Summer School in Glaciology](http://glaciers.gi.alaska.edu/courses/summerschool), McCarthy, AK in years 2010, 2012, 2014, 2016, 2018, and 2022.  The PDF notes (the `.pdf` in [notes/](notes/)) plus their exercises plus the Matlab/Octave programs form a self-contained course.  The PDF slides (the `.pdf` in [slides/](slides/)) cover the same material in a more informal style.


download Matlab/Octave programs
-------------------------------

The codes in subdirectory [mfiles/](mfiles/) solve SIA and SSA problems.  They should work in either Matlab or [Octave](https://www.gnu.org/software/octave/); if not please report a bug, either by email or by using the [issues](https://github.com/bueler/mccarthy/issues) for this repository.  Download them either by cloning this repository or by getting a `.zip` or `.tar.gz` archive at the [releases page](https://github.com/bueler/mccarthy/releases).  The PDFs mentioned above document these programs, but the programs also have help files (i.e. leading comments).  You are encouraged to actually run and modify them!


stokes solver
-------------

The Python tools in [stokes/](stokes/) solve a free-surface 2D Glen-Stokes flow over a bedrock step.  This is used to generate some images in the slides and for student projects.  Note that the numerical technology used here is more advanced, and more experimental, than the [mfiles/](mfiles/) content.  The workflow uses the following tools: [Firedrake finite element library](https://www.firedrakeproject.org/), [Gmsh mesh generator](http://gmsh.info/), [PETSc solver library](http://www.mcs.anl.gov/petsc/), and [Paraview visualization](https://www.paraview.org/).  See [stokes/README.md](stokes/README.md) for more information.


versions
--------

Older versions (2009, 2010, 2012, 2014) lived in the repo https://github.com/bueler/karthaus
