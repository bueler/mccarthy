mccarthy
========

Copyright 2010--2026  Ed Bueler

This repository contains slides, notes, exercises, and computer programs for numerical glacier and ice sheet modeling.  These materials have been used for the [International Summer School in Glaciology](http://glaciers.gi.alaska.edu/courses/summerschool), McCarthy, AK in years 2010, 2012, 2014, 2016, 2018, 2022, 2024, and 2026.

Download
--------

Download all materials by one of these methods:

  * use the [releases page](https://github.com/bueler/mccarthy/releases) for a `.zip`/`.tar.gz` archive
  * clone or shallow clone this repository:

    git clone --depth=1 https://github.com/bueler/mccarthy.git

The optional `depth=1` setting reduces the download size by not getting the history, which most users will not need.

Slides, exercises, and notes
----------------------------

The PDF slides [slides/slides-2026.pdf](slides/slides-2026.pdf) cover the essential material in the briefest way.  They are supported by Python codes in [py/](py/).  There are [exercises](exercises/exercises-2026.pdf) and [notes](notes/notes-2026.pdf) which are aligned to the slides, and also support them.

Python programs
---------------

The Python codes in subdirectory [py/](py/) solve surface kinematic equation (SKE) and shallow ice approximation (SIA) problems.  The primary Python programs [surface1d.py](py/surface1d.py) and [shallowuw.py](py/shallowuw.py) use only the numpy and matplotlib libraries.

You are encouraged to actually run and modify these codes!  Feel free to email me about how they work.

The [slides](slides/slides-2026.pdf) and the [notes](notes/notes-2026.pdf) are the documentation for the Python codes, along with line comments in the `.py` source codes.

The [surface2d.py](py/surface2d.py) program uses [Firedrake](https://www.firedrakeproject.org/) (a finite element library) which calls [PETSc](http://www.mcs.anl.gov/petsc/) for equation solvers.

The codes in [py/stokes/](py/stokes/) also use Firedrake to solve a 2D Glen-law Stokes flow over a bedrock step.  The workflow for the Stokes codes also use [Gmsh](http://gmsh.info/) (a mesh generator) and [Paraview](https://www.paraview.org/) (for visualization).  See [py/stokes/doc.pdf](py/stokes/doc.pdf) for more information.

Please report any bugs in the codes either by email or by using the [issues](https://github.com/bueler/mccarthy/issues) for this repository.

Older material
--------------

Older and deprecated material is in the form of PDF notes (see [old/notes/](old/notes/)), slides (see [slides/old/](slides/old/), and Matlab/Octave programs [old/mfiles/](old/mfiles/).  Even older versions (2009--2014) of this material lives in the repo [github.com/bueler/karthaus](https://github.com/bueler/karthaus).
