Glacier flow over a cliff
=========================

reassurance
-----------

The point of the project is not to reach a particular destination!  The point is to _help you learn fast from whereever you currently are_.  So please be honest with yourself about what you don't know, and ask about it.  Another student may know more programming than you, or have seen more differential equations, or have more glaciers background, or whatever.  But things like that are irrelevant to my goal, which is to help you move forward as far as possible in one intense week.  Forward progress must be from the things _you_ understand into things _you_ do not understand.


description
-----------

PROJECT 4: Glacier flow over a cliff

ADVISOR: Ed Bueler

DESCRIPTION: How does a glacier flow over a submerged cliff or step in the bedrock?  This geometry can be put into in a complete Stokes flow model.  We will work with one of such a model and thereby get a pretty good answer.  On the other hand, the shallow models used for large-scale ice sheet studies must also produce reasonable results for bedrock configurations including such steep cliffs and other bad bedrock topography.  This project will work with the equations and numerical models associated to both Stokes and shallow stress balances.  We aim to build an understanding of the flow and glacial profiles associated to bad topography, and perhaps figure out good ways to modify a shallow model so as to handle bed cliffs effectively.

SOFTWARE REQUIREMENTS: Numerical Python (Python+numpy+matplotlib).  The Firedrake finite element Python library will also be used, but that will come installed on a provided laptop.

REQUIRED STUDENT BACKGROUND: Exposure to differential equations and linear algebra.  Basic programming with numerical Python.


references and reading
----------------------

See the notes `notes/notes-bueler-2018.pdf` in the current repository.  This project directly extends the material in section 2 regarding the planar Stokes model and the slab-on-slope case.  Equally relevant is the derivation of the SIA (shallow ice approximation) model in sections 2 and 3 of the notes.

A closely-related study of the differences between planar Stokes flow and the SIA appears in

  * G. J. Leysinger Vieli and G. H. Gudmundsson, 2004.  _On estimating length fluctuations of glaciers caused by changes in climatic forcing_, J. Geophys. Res.: Earth Surface 109, F01007, https://doi.org/10.1029/2003JF000027.

This article is available in the current directory (and possibly on paper).  Leysinger-Vieli and Gudmundson show that the SIA is actually rather good at tracking changes in ice geometry, even with far-from-shallow ice geometry, assuming smooth bed topography.  I'm curious what is the best way to fix the SIA to get results which are close to the Stokes solution.

