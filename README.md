mccarthy
========

Copyright 2010--2018  Ed Bueler

These are slides, notes, and codes on numerical glacier and ice sheet modeling, for the International Summer School in Glaciology, McCarthy, AK.  The notes (`notes/notes.pdf`) plus their exercises plus the Matlab/Octave codes form a self-contained course.  The slides (`slides/lecture.pdf`) cover the same material in a more informal style.

The Summer School website is at http://glaciers.gi.alaska.edu/courses/summerschool


download Matlab/Octave codes
----------------------------

The codes in subdirectory `mfiles/` should work both in Matlab and [Octave](https://www.gnu.org/software/octave/); if not please report a bug, by email or using the [issues](https://github.com/bueler/mccarthy/issues) for this repo.  Download these codes either by cloning this repo or by getting a "release" in `.zip` or `.tar.gz` archive format at https://github.com/bueler/mccarthy/releases
and unpacking it.

The PDFs above (`notes/notes.pdf` and `slides/lecture.pdf`) are the major documentation of these programs, but the codes also have help files (i.e. leading comments).  You are encouraged to actually run and modify the codes!


stokes solver
-------------

The Python tools in `stokes/` are used to generate some images in the slides, and for some student projects.  Such numerical technology is well outside the scope of the notes and slides; it is more advanced and more experimental.  They exploit the Firedrake finite element framework (https://www.firedrakeproject.org/), the Gmsh mesh generator (http://gmsh.info/), the PETSc solver library (http://www.mcs.anl.gov/petsc/).  See `stokes/README.md` for more information.


versions
--------

This 2018 version updates the 2016 version.  Older versions (2010, 2012, 2014) lived in the repo https://github.com/bueler/karthaus

