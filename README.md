mccarthy
========

Copyright 2010--2018  Ed Bueler

These are slides, notes, and codes on numerical glacier and ice sheet modeling, for the International Summer School in Glaciology, McCarthy, AK.  The notes (`notes/notes.pdf`) plus their exercises plus the Matlab/Octave codes form a self-contained course.  The slides (`slides/lecture.pdf`) cover the same material in a more informal style.

The Summer School website is at

http://glaciers.gi.alaska.edu/courses/summerschool


download Matlab/Octave codes
----------------------------

The codes are in subdirectory `mfiles/`.  Download these codes either by cloning this repo or by getting a "release" in `.zip` or `.tar.gz` archive format at

        https://github.com/bueler/mccarthy/releases

and unpacking it.

The PDFs above are the major documentation of these programs, but they also have help files (i.e. leading comments).  You are encouraged to actually run and modify the codes!


stokes solver
-------------

The Python tools in `stokes/` are used to generate some images in the slides, and for some student projects.  They exploit the Firedrake finite element framework (https://www.firedrakeproject.org/), the Gmsh mesh generator (http://gmsh.info/), and the PETSc solver library (http://www.mcs.anl.gov/petsc/).  Such numerical technology is well outside the scope of my notes and slides for the Summer School.  Do not expect good introductory documentation when looking in this directory.


versions
--------

This is the 2018 version, which is an update of the 2016 version.  Older versions (2010, 2012, 2014) lived inside the repo https://github.com/bueler/karthaus

