# Cliffs, overhangs, damage

## reassurance

The point of the project is not to reach a particular destination!  It is to _help you learn fast from wherever you currently are_.  So please be honest with yourself about what you don't know, and please ask about it!  Another student may know more programming than you, or have seen more differential equations, or have more glaciers background, or data-manipulation skills, or whatever.  That is irrelevant to our goals.  We just want to help you move forward as far as possible in one intense week.  Forward progress must be from things _you_ understand into things _you_ do not yet understand!

## description

PROJECT 10: Cliffs, overhangs, damage

ADVISOR: Ed Bueler

DESCRIPTION: **FIXME** How does the velocity field of a glacier change if you put a pile of ice somewhere on its surface, or if you dig a pit somewhere?  The change in the velocity field reveals the boundary Green's functions of the Stokes problem for glaciers.  The goal of this project is to build new understanding of these Green's functions.  How do they depend on glacier thickness or sliding?  How do they depend on nearness to margins or bed topography?  How can they be exploited to build faster numerical models?  Students will use a modern finite-element framework for solving the Glen-Stokes equations, in a flowline plane, for student-chosen glacier geometries, as we explore.

SOFTWARE REQUIREMENTS: Python _running locally on your machine_.  Please try to build/install: Firedrake, Gmsh, Paraview.  (As backup, I'll bring an extra laptop, pre-loaded.)

REQUIRED STUDENT BACKGROUND: Exposure to linear algebra, numerical methods, partial differential equations, and (perhaps) the finite element method.

## getting started

As the first step, clone the current [McCarthy repo](https://github.com/bueler/mccarthy):

    $ git clone --depth=1 https://github.com/bueler/mccarthy.git

Go to the `stokes/` directory and try all the steps documented in `stokes/README.md`.  You will learn how to use my Firedrake-based Python Stokes solver, with pre- and post-processing using Gmsh and Paraview.  Please also get started on reading the documentation `stokes/stokes.pdf` for the solver.  Some of these steps may be quite new to you, in which case the project has already been worthwhile!

My suggested goal of this project is to **FIXME**