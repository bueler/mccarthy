## PROJECT 8: How far can glacier stresses reach?

ADVISOR: Ed Bueler

DESCRIPTION:  How far does a localized bed topography perturbation, or a localized perturbation in basal resistance, influence the glacier velocity?  What controls the range of that influence?  Students will explore which physical parameters (e.g. aspect ratio, slip ratio, rheology, temperature profile) extend or shrink the range of influence in a Stokes model.  The finite element method will be used to solve a flow-line Stokes stress balance model, and students will learn about meshing and solver choices as well.

Required tools:  This project assumes some familiarity with Python and differential equations, a knowledge basis which we will strengthen.  Students will learn the Firedrake finite element library, the Gmsh meshing program, and the Paraview tool for visualization.

### Getting started

Please read `py/stokes/doc.pdf` lightly.  (I brought printed-out versions for low-tech McCarthy usage, and I'll hand them out.)  Then run the example at the end.  This will make sure that you have seen a geometry outline (`.geo` file), how Gmsh generates a mesh, how the main Firedrake program (`py/flow.py`) runs and solves a Stokes problem, and then have seen how Paraview visualizes the solution.  Please debug any issues at this stage, and seek my help with anything that comes up.  I hope you will talk to each other as you do this.

Let's talk, as much as you can stand, about the ideas in these codes.  We can talk through how the planar Stokes equations work; this is also in my regular notes.  I am happy to discuss, as much as you want, the process of solving partial differential equations by the finite element method, via the weak form and Newton's method; these steps are not obvious historically-speaking, and they take everyone some time to learn.

For the project itself, I will try to help you _explore the ideas which you find most interesting_.  My specific goal here is to communicate an idea: the glaciological stress balance equations can be understood through making small, local perturbations at points.  These perturbations can be implemented through small changes to my codes that you can make; they only need a few lines changed.  To see the effect of a perturbation first solve for the background state.  Then change a value at one node, rerun, and subtract the background field from the new (perturbed) field.  Visualize and measure the difference.

The question of how far stresses reach from such a perturbation, in the original project description above, is just one version of what numerical perturbations can reveal.  You are not limited to doing that; feel free to pursue what you find interesting.

### Formulating a plan

You will want to focus on some particular kind of perturbation.  There are certainly others one could consider, but here are some options for that:

1. Make a small change to surface topography or base topography by slightly-modifying the mesh coordinates of one node.

2. Make a small change to the upstream input velocity, or to the downstream stress condition.

3. Implement a patch of stress-free base, and make small perturbations to its extent.

4. Implement a patch of sliding-law base, and make small perturbations to its extent.

5. The original project description mentions changes to: aspect ratio, slip ratio, rheology, temperature profile.  (The last two would presumably be through parameters in the flow law.)

In any case, modifying and experimenting with the codes, to make them do new things, is a route to understanding their underlying ideas.  For the project presentation, please try to formulate a hypothesis which your modifications and experiments will test.

### Reassurance

It may turn out that a lot of what you do comes down to learning some skills which are necessary background for numerical glacier modeling.  Such as: comfort editing and running programs from the command line, learning how to program in Python, learning how finite difference and/or finite element numerical methods work, learning how meshing works, or learning how Firedrake works.  In which case, this is good!  These are useful skills!
