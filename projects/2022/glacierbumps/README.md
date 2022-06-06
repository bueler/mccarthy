# Stokes modeling of glacier velocity: the effects of surface bumps

## reassurance

The point of the project is not to reach a particular destination!  It is to _help you learn fast from wherever you currently are_.  So please be honest with yourself about what you don't know, and please do ask about it!  Another student may know more programming than you, or have seen more differential equations, or have more glaciers background, or data-manipulation skills, or whatever.  That is irrelevant to my goal.  I just want to help you move forward as far as possible in one intense week.  Forward progress must be from things _you_ do understand into things _you_ do not yet understand!

## description

PROJECT 11: Stokes modeling of glacier velocity: the effects of surface bumps

ADVISOR: Ed Bueler

DESCRIPTION: How does the velocity field of a glacier change if you put a pile of ice somewhere on its surface, or if you dig a pit somewhere?  The change in the velocity field reveals the boundary Green's functions of the Stokes problem for glaciers.  The goal of this project is to build new understanding of these Green's functions.  How do they depend on glacier thickness or sliding?  How do they depend on nearness to margins or bed topography?  How can they be exploited to build faster numerical models?  Students will use a modern finite-element framework for solving the Glen-Stokes equations, in a flowline plane, for student-chosen glacier geometries, as we explore.

SOFTWARE REQUIREMENTS: Python/numpy, _running locally on your machine_.  Ideally, functional Firedrake, Gmsh, and Paraview installations.  (But I'll bring an extra laptop, pre-loaded.)

REQUIRED STUDENT BACKGROUND: Exposure to linear algebra, numerical methods, partial differential equations, and (ideally) the finite element method.

## getting started

As the zeroth step, clone [this McCarthy repo](https://github.com/bueler/mccarthy) and go to the `stokes/` directory and make sure you can run all the steps documented in `stokes/README.md`, which describe how to use my Firedrake-based Python Stokes solver, with pre/post-processing using Gmsh and Paraview.  (This may already be a project, which is fine!)  Also get started on reading the documentation `stokes/stokes.pdf` for the solver.

My suggested goal of this project is to compare, *visually and quantitatively*, three different meanings of the *glaciological Green's functions for surface changes*.  These Green's functions, also called _Stokeslets_, are solutions of a linear problem, namely the linearization of the full glaciological Glen-Stokes problem for a choice of glacier geometry and boundary stresses.  We can compute them numerically by using my existing Glen-Stokes finite element solver and applying it for slightly-perturbed problems.

Remember that our only access to Glen-Stokes solutions, for any other geometry than slab geometry, is through numerics.  Thus my explanation here is in terms of a finite element mesh and solutions thereon.  For simplicity I'll only consider the frozen-base case, with zero Dirichlet boundary condition at the base ($\mathbf{u}=0$), and I'll assume a stress-free condition at the surface (so $p=0$ there).

Suppose we set up some glacier geometry and compute the (numerical) solution to the Glen-Stokes equations.  This solution is not the Green's function; this is the _background field_ $\mathbf{u}$ which we now perturb.

Here are three ways to explain physically what are the Green's function for a given location $x$ on the surface of a glacier.  Each meaning describes a change in the solution relative to the background field:

  1. The _change_ in the background field caused by adding a unit mass of ice (e.g.~one tonne, $10^3$ kg) to the surface at $x$ by moving a mesh node location upward by an appropriate amount so that the added volume corresponds to the unit mass.  (The mesh will change but the equations and boundary conditions will remain the same.)
  2. The _change_ in the background field caused by increasing the force of gravity "at" the mesh node at location $x$ by an amount equivalent to the. One adds this mass by increasing the mass of the elements which meet that node.  (The right-hand side of the equations will change but the mesh and the boundary conditions will remain the same.)
  3. The _change_ in the background field caused by adding a nonzero boundary stress equivalent to the weight of a unit mass applied at location $x$.  (The boundary conditions will change but the mesh and the right-hand side will remain the same.)

After discussing strategies for actually generating these fields, I will encourage you to choose your own route to computing and visualizing these fields.

There are at least two background geometries which are set to go.  They correspond to Python scripts which generate `.geo` files for mesh generation by Gmsh:

  * `stokes/domain.py` generates a geometry which is a slab with a bedrock bump
  * `stokes/project11/profile.py` generates a geometry for an ice sheet margin

I am imagining the generation of a visual gallery of Green's functions, as well as a more quantitative study of the Green's functions as they act back at the surface.  There is much to discuss, and many directions to go, and you will be very much in charge of which aspects you want to pursue.
