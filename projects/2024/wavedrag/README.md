# Drag from waves in the glacier bed

## reassurance

The point of the project is not to reach a particular destination.  It is to _help you learn fast from wherever you currently are_.  So please be honest with yourself about what you don't know, and please ask about it!  Another student may know more programming than you, or have seen more differential equations, or have more glaciers background, or better data manipulation skills, or whatever.  That is irrelevant to our goals.  We want to help you move forward as far as possible in one intense week.  Forward progress will be from things _you_ understand into things _you_ do not yet understand.

## description

PROJECT 13: Drag from waves in the glacier bed

ADVISOR: Ed Bueler

DESCRIPTION:  How much drag is caused by waves in the bed topography of a glacier?  How does this depend on wavelength and magnitude?  Does sliding make the effect larger or smaller?  These questions can be explored using a numerical Stokes model.  We will compare slab-on-slope solutions to cases where the base has added topography or sliding, and measure the modeled velocity change.  Because these questions have been addressed by analytical/expansion techniques in earlier literature, we can also use the model results to assess those conclusions.  Students will gain experience and understanding as they modify an already-written finite-element solver of the 2D (planar) Glen-Stokes equations.

SOFTWARE REQUIREMENTS: You will need a recent version of Python _running locally on your machine_.  Please try to build/install the following: Firedrake, Gmsh, Paraview.  (As backup, I'll bring an extra laptop, pre-loaded.)

STUDENT BACKGROUND: Linear algebra and some differential equations are required.  Optionally, perhaps a bit of numerical methods or the finite element method?

## getting started

As the first step, if you have not done it already, please clone my whole [McCarthy repo](https://github.com/bueler/mccarthy):

    $ git clone --depth=1 https://github.com/bueler/mccarthy.git

Now go to the `stokes/` directory and try all the steps documented in `stokes/README.md`.  You will learn how to use my Firedrake-based Python Stokes solver, with pre- and post-processing using Gmsh and Paraview.  Some of these "first step" steps may involve some learning, in which case the project has already been worthwhile!

As a second step, collect and browse the following documents.  I will provide these on paper in McCarthy, but they are also available electronically:

  * `stokes/doc.pdf`:  This documents how the Firedrake Stokes solver works.  It can be built in `doc/` using LaTeX.
  * `projects/2024/wavedrag/Schoof2003`:  This paper contains FIXME
