# projects/2026/

These are notes and codes for my two suggested projects for the McCarthy International Summer School in Glaciology in 2026.  See `README.md` in the subdirectories `coverage/` (Project XX) and `influence/` (Project YY) for documentation.  These projects appear in the [projects list](https://glacierschool.alaska.edu/2026/student-projects-2026) for the Summer School.

## coverage/

PROJECT XX: What land is covered by glaciers?

ADVISOR: Ed Bueler

DESCRIPTION:  Given a _prescribed_ ice surface velocity field, a general distribution of surface mass balance, and a more-or-less realistic glacier bed, what land is covered by glaciers?  The severely-simplified model behind this question can be steady-state or time-dependent, but it is a free-boundary problem in which the ice thickness is unknown.  Three fields, namely surface mass balance, bed topography, and the (strangely) predetermined velocity field, determine glacier coverage and thickness (geometry).  We will solve this model numerically using a finite element model.  Then variations will be entertained.  Students might compare two formulations, thickness-based mass continuity and the surface kinematic equation.  Students might improve the velocity field toward something physically-believable, such as a shallow approximation of the glacier flow physics.

Required tools:  This project assumes some familiarity with Python and differential equations, a knowledge basis which we will strengthen.  Students will learn the Firedrake finite element library and the Paraview tool for visualization.

## influence/

PROJECT YY: How far can glacier stresses reach?

ADVISOR: Ed Bueler

DESCRIPTION:  How far does a localized bed topography perturbation, or a localized perturbation in basal resistance, influence the glacier velocity?  What controls the range of that influence?  Students will explore which physical parameters (e.g. aspect ratio, slip ratio, rheology, temperature profile) extend or shrink the range of influence in a Stokes model.  The finite element method will be used to solve a flow-line Stokes stress balance model, and students will learn about meshing and solver choices as well.

Required tools:  This project assumes some familiarity with Python and differential equations, a knowledge basis which we will strengthen.  Students will learn the Firedrake finite element library, the Gmsh meshing program, and the Paraview tool for visualization.
