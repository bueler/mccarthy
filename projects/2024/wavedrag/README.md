# FIXME

## reassurance

The point of the project is not to reach a particular destination!  It is to _help you learn fast from wherever you currently are_.  So please be honest with yourself about what you don't know, and please ask about it!  Another student may know more programming than you, or have seen more differential equations, or have more glaciers background, or data-manipulation skills, or whatever.  That is irrelevant to our goals.  We just want to help you move forward as far as possible in one intense week.  Forward progress must be from things _you_ understand into things _you_ do not yet understand!

## description

PROJECT 11: FIXME

ADVISOR: Ed Bueler

DESCRIPTION:  FIXME At the surface of a glacier, and especially on steep margins, our viscous-fluid understanding of glaciers can break down.  Near cliffs and overhangs, stresses within the ice may turn into fractures and crevasses.  A numerical Stokes model can address how fractures appear via a model of damage, that is, of the deterioration of polycrystalline structure.  We will try to model the initial damage, starting from a good base in the literature and an already-written, and brief, finite-element solver of the 2D (planar) Glen-Stokes equations.  Applying the solver will connect ice geometry and surface stresses to stresses within the ice.  These stresses can then be expressed as rates of change of damage.  We will explore different geometries, evolution models, and questions as they arise.

SOFTWARE REQUIREMENTS: You will need a recent version of Python _running locally on your machine_.  Please try to build/install the following: Firedrake, Gmsh, Paraview.  (As backup, I'll bring an extra laptop, pre-loaded.)

STUDENT BACKGROUND: Linear algebra and some differential equations are required.  Optionally, perhaps a bit of numerical methods or the finite element method?

## getting started

As the first step, if you have not done it already, please clone my whole [McCarthy repo](https://github.com/bueler/mccarthy):

    $ git clone --depth=1 https://github.com/bueler/mccarthy.git

Now go to the `stokes/` directory and try all the steps documented in `stokes/README.md`.  You will learn how to use my Firedrake-based Python Stokes solver, with pre- and post-processing using Gmsh and Paraview.  Some of these "first step" steps may involve some learning, in which case the project has already been worthwhile!

FIXME As a second step, collect and browse the following documents.  I will provide these on paper in McCarthy, but they are also available electronically:

  * `stokes/doc/stokes.pdf`:  This can be built using LaTeX.  It documents how the Firedrake Stokes solver works.
  * `projects/2024/wavedrag/Schoof2003`:  This paper contains FIXME

My suggested goal of this project is to **FIXME**