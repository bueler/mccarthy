## PROJECT 7: What land is covered by glaciers?

ADVISOR: Ed Bueler

DESCRIPTION:  Given a _prescribed_ ice surface velocity field, a general distribution of surface mass balance, and a more-or-less realistic glacier bed, what land is covered by glaciers?  The severely-simplified model behind this question can be steady-state or time-dependent, but it is a free-boundary problem in which the ice thickness is unknown.  Three fields, namely surface mass balance, bed topography, and the (strangely) predetermined velocity field, determine glacier coverage and thickness (geometry).  We will solve this model numerically using a finite element model.  Then variations will be entertained.  Students might compare two formulations, thickness-based mass continuity and the surface kinematic equation.  Students might improve the velocity field toward something physically-believable, such as a shallow approximation of the glacier flow physics.

Required tools:  This project assumes some familiarity with Python and differential equations, a knowledge basis which we will strengthen.  Students will learn the Firedrake finite element library and the Paraview tool for visualization.

### Getting started, and a formulating a plan

1. Read and run `py/surface.py`, which is in 1D, and read and run `py/surface2d.py` also.  The more elementary numerical methods in `surface.py`, which does _not_ use Firedrake, are emphasized in my slides and notes.  (It is possible to do a nice project based on this 1D code, and not use Firedrake at all.)  The ideas in `surface2d.py`, for 2 horizontal dimensions (2D), include several advanced concepts, including weak forms, variational inequalities, and solvers.  Talk to each other as you do this.

2. Let's talk, as much as you can stand, about the ideas in these codes.  I will try to help you _explore the ideas which you find most interesting_.  My goal for this project is to communicate the idea that glaciers end where ablation has driven their thicknesses to zero, as they flow downhill.  This is an interesting free-boundary problem; understanding how it works is part of designing good numerical glacier models.

3. The basic outline of the project is to modify and experiment with the codes, and make them do new things, as a route to understanding their underlying ideas.  Please try to formulate a conceptual hypothesis which your modifications and experiments will test.  Options for how to proceed include:

    a) Add ice dynamics to the 1D code `surface.py` in its steady-state mode.  It could be converted into an shallow ice approximation (SIA) glacier simulation.  This could be done by inventing an iteration scheme, alternating between solving the surface kinematical equation (SKE).  In such an iteration you would apply `surface.py` to convert velocity into glacier geometry, as it currently does, and then compute the velocities from the SIA and the geometry, which is what `shallow.py` currently does.  Then you would try to make this iteration converge.

    b) Modify the 2D code `surface2d.py`, which is also steady state, to have more interesting, or at least significantly different, data.  The data are: topography, surface mass balance, and velocity.  The data could be randomly modified, or there could be specific features to change.  Numerical experiments could study where the free boundary goes under such perturbations.

    c) Modify the 2D code `surface2d.py` to have ice dynamics from the SIA.  The easiest modification is to do this in a time-dependent manner.  This is a coding-intensive path.

    d) Modify the 1D code `surface.py` to have ice dynamics from the Stokes equations.  Again, the easiest modification is to do this in a time-dependent manner.  This is also a coding-intensive path; it would involve learning what project 8 is doing.

4. It may turn out that a lot of what you do comes down to learning some skills which are necessary background for numerical glacier modeling.  Such as: comfort editing and running programs from the command line, learning how to program in Python, learning how finite difference and/or finite element numerical methods work, or learning how Firedrake works.  In which case, this is good!  These are useful skills!
