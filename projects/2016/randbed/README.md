
Glaciating random terrain
=========================


description
-----------

The original description:

> PROJECT 8: Glaciating random terrain

> STUDENTS: XXX

> ADVISOR: Ed Bueler

> DESCRIPTION: How does glacier-covered area, or glacier volume, depend on the roughness of topography?  In this modeling project we will generate random terrain, assume altitude-accumulation feedback, model the corresponding steady states of the glaciers, and measure the modeled glaciated area and glacier volume.  How does glaciation change as equilibrium line altitude and/or lapse rates change on random terrain of different steepness?  This project will involve running laptop-suitable glacier models and processing their output.

> SOFTWARE REQUIREMENTS: Python and Linux/Unix

> REQUIRED STUDENT BACKGROUND: Basic exposure to (1) differential equations, (2) some numerical methods, and (3) use of Python and unix-type command line.


references
----------

First, see the notes (`notes/notes.pdf`).  Pay special attention to the parts about the 2D SIA model.  For that model there are issues to discuss about degenerate diffusivity and the fact that the thickness is constrained to be nonnegative.

Second, this is available on paper and in the current directory: E. Bueler, 2016.  _Stable finite volume element schemes for the shallow-ice approximation_, J. Glaciol. (published online) doi:10.1017/jog.2015.3.

Third, regarding generating a random bed, see https://en.wikipedia.org/wiki/Gaussian_process.


needed tools
------------

_Ideally_ all tools would work on your laptop, but we can manage if not.  The first two PETSc items may be hard to install.  The tools are:

* PETSc.  See http://www.mcs.anl.gov/petsc/ and follow download and install instructions.

* My SIA-solving program, written in C and using PETSc, in support of Bueler (2016) above.  Clone my repo https://github.com/bueler/sia-fve.  The code is documented by Bueler (2016) above, except for the elevation-dependence of the climatic mass balance (CMB).  This part is  m(x,y) = lapse * (s(x,y) - ela)  where m(x,y) is the CMB (m s^-1) and where  s(x,y) = b(x,y) + H(x,y)  is the surface elevation (m).

* python with numerical and graphical libraries.  See https://www.python.org/, http://www.scipy.org/, http://matplotlib.org/.

* My random-bed-generation program in the current directory: `bedrand.py`.  Looking at the 1D Matlab version `bedrand1d.m` may help understand the basics.

* Any tools you need to post-process results and produce a short talk at the end of the school.


questions
---------

* What is the continuum model?  Be able to state equations and explain what they mean.

* What are the differences between the time-dependent and steady state continuum models?

* What kinds of things are known and not known about the model?  (Solutions?  Numerical solutions?  Appropriateness of the continuum model?)

* Someone might say that the steady state may not be unique if the CMB is elevation-dependent.  What does this mean?

* How do you measure the extent/magnitude/coverage of the glaciation of a mountain range?

* What are the parameters controlling the random bedrock topography?  What cases/values should be studied?

* What are the parameters controlling the continuum ice flow model?  What cases/values should be studied?

* What are the additional parameters controlling the _numerical_ ice flow model?  What cases/values should be studied?

* Why do lower resolution runs converge less reliably than higher resoluion?  (See below for examples of lower and higher.  This is a hard question ... I do not really know but I can guess!)


goals/mileposts
---------------

* Decide on an initial experimental design _quickly_.  That is, understand the basic idea and the basic models, and then take a guess at what to look at.

* Key idea: _I can help run the experiments._  I can be a flunky lab assistant who knows how to run the computer programs and wait for them to finish and collect the outputs on e.g. a USB stick.

* Generate the first part of your presentation early.  I.e. slides stating the model, the parameters, the goals, and the initial experiment.  (You have to do this anyway, and doing it early helps reveal good questions early; that is the goal!)


recipes
-------

Build PETSc code

        (cd ~/sia-fve/petsc/ && make mahaffy)

Generate links:

        ln -s ~/sia-fve/petsc/mahaffy
        ln -s ~/sia-fve/petsc/figsmahaffy.py
        ln -s ~/petsc/bin/PetscBinaryIO.py
        ln -s ~/petsc/bin/petsc_conf.py

Generate and view example bed:

        ./bedrand.py --plotbed

Generate/view, and then save, example bed:

        ./bedrand.py --plotbed -o low.dat

Run a low-resolution example using the bed:

        mkdir testlow
        ./mahaffy -mah_read low.dat -mah_showdata -draw_pause 2 -mah_cmbmodel -cmb_ela 2500.0 -cmb_lapse 0.002 -cs_D0 1.0 -snes_max_it 400 -mah_dump testlow/

Look at result:

        cd testlow/
        ../figsmahaffy.py --observed             # view resulting .png and .pdf files

Run a higher-resolution example, in parallel:

        ./bedrand.py --plotbed -Nx 100 -Ny 100 -o high.dat
        mkdir testhigh
        mpiexec -n 4 ./mahaffy -mah_read high.dat -mah_showdata -draw_pause 2 -mah_cmbmodel -cmb_ela 2500.0 -cmb_lapse 0.002 -cs_D0 1.0 -mah_dtrecovery 10 -snes_max_it 400 -pc_type asm -sub_pc_type lu -mah_dump testhigh/
        cd testhigh/
        ../figsmahaffy.py --observed

