## PROJECT 8: How far can glacier stresses reach?

ADVISOR: Ed Bueler

DESCRIPTION:  How far does a localized bed topography perturbation, or a localized perturbation in basal resistance, influence the glacier velocity?  What controls the range of that influence?  Students will explore which physical parameters (e.g. aspect ratio, slip ratio, rheology, temperature profile) extend or shrink the range of influence in a Stokes model.  The finite element method will be used to solve a flow-line Stokes stress balance model, and students will learn about meshing and solver choices as well.

Required tools:  This project assumes some familiarity with Python and differential equations, a knowledge basis which we will strengthen.  Students will learn the Firedrake finite element library, the Gmsh meshing program, and the Paraview tool for visualization.

### Getting started, and a formulating a plan

1. I have printed out `py/stokes/doc.pdf`.  Please read it lightly.  Then run the example at the end.  This will make sure that you have seen Gmsh generate a mesh, then have the main Firedrake program (`py/flow.py`) solve a Stokes problem, and then use Paraview to view the answer.  Debug any issues at this stage, and seek my help with any stage.  Talk to each other as you do this.

2. Let's talk, as much as you can stand, about the ideas in these codes.  In particular, we can talk through the planar Stokes equations, which are also in my regular notes.  Also the process of solving partial differential equations by the finite element method, via the weak form and Newton's method, is not obvious.  I want to help you understand!  For the project itself, I will try to help you _explore the ideas which you find most interesting_.

3. My specific goal here is to communicate an idea: the glaciological Stokes equations can be understood through making small, local perturbations at points.  These perturbations can be implemented through small changes to my codes that you can make; they only need a few lines changed.  The key is to solve once for the background state, then change a value at one node, then rerun, then subtract the background field from the changed (perturbed) field, and look at the difference.  The question of how far stresses reach, from such a perturbation, is just one version of what numerical perturbations can reveal.

4. You will want to focus on some particular kind of perturbation.  There are certainly others one could consider, but here are some options for that:

    a) Make a small change to surface topography or base topography by slightly-modifying the coordinates of one point.

    b) Make a small change to the upstream input velocity, or to the downstream stress condition.

    c) Implement a patch of stress-free base, and make small perturbations to its extent.

    d) Implement a patch of sliding-law base, and make small perturbations to its extent.

    e) The original project description at the top mentions changes to aspect ratio, slip ratio, rheology, temperature profile.

5. Modifying and experimenting with the codes, and to make them do new things, is a route to understanding their underlying ideas.  However, for the project presentation, please try to formulate a hypothesis which your modifications and experiments will test.

6. It may turn out that a lot of what you do comes down to learning some skills which are necessary background for numerical glacier modeling.  Such as: comfort editing and running programs from the command line, learning how to program in Python, learning how finite difference and/or finite element numerical methods work, learning how meshing works, or learning how Firedrake works.  In which case, this is good!  These are useful skills!
