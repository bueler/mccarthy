# projects/2026/

These are notes and codes for my two suggested projects for the McCarthy International Summer School in Glaciology in 2026.  See `README.md` in the subdirectories `coverage/` (Project XX) and `influence/` (Project YY) for documentation.  These projects appear in the [projects list](https://glacierschool.alaska.edu/2026/student-projects-2026) for the Summer School.

## coverage/

Given a _prescribed_ ice surface velocity field and a more-or-less realistic glacier bed, determine where glaciation covers the land.  This severely-simplified problem can be steady-state or time-dependent, but the point is that it is posed as a free-boundary problem in which the ice thickness is unknown and the surface mass balance and a velocity field determine glacier coverage (and thickness).  Students might compare two formulations, thickness-based mass continuity and the surface kinematic equation.  Students may refine the velocity field toward something physically-believable, such as a shallow approximation of the glacier flow physics.

Required tools:  This project assumes some familiarity with Python and differential equations, which is a basis we will strengthen.  Students will learn the Firedrake finite element library and the Paraview tool for visualization.

## influence/

FIXME How far does a localized bed perturbation influence the glacier surface, and what controls the range of that influence?  This is related to `projects/2022/glacierbumps/` but the focus is not on the Green's function interpretation.  Instead, students explore what physical parameters (e.g. slip ratio, geometry, rheology) extend or shrink the range of influence in a Stokes model.

Required tools:  This project assumes some familiarity with Python and differential equations, which is a basis we will strengthen.  Students will learn the Gmsh meshing program, the Firedrake finite element library and the Paraview tool for visualization.
