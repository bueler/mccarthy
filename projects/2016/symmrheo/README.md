
Ice rheology from symmetric flows
---------------------------------


description
-----------

The original description:

> PROJECT 7:  Ice rheology from symmetric flows

> STUDENTS: XXX

> ADVISOR: Ed Bueler

> DESCRIPTION:  There are simple models for the evolution of circularly-symmetric flows of ice (or other viscous fluids).  This theory, used quantitatively, allows us to convert observations of real flows, whether approximately circular (polar caps on Mars) or deliberately circular (laboratory), into estimates of the rheology of the fluid.  The project will use pencil and paper calculations, some processing of data on computer, and modest experimentation.

> SOFTWARE REQUIREMENTS: Matlab or Octave

> REQUIRED STUDENT BACKGROUND: Basic exposure to (1) differential equations and (2) Matlab or similar.


references
----------

First, see the notes (`notes/notes.pdf`).  They cover the shallow ice approximation (SIA; equation (26) in notes) which you will want to comprehend, and the Halfar (1983) similarity solution (equation (40) in notes).

The Halfar solution is used in the following paper, a PDF of which is in the current directory: J. F. Nye, 2000. _A flow model for the polar caps of Mars_, J. Glaciol. 46 (154), 438--444.

The second paper to read is the one for which I am proposing you reproduce part of the analysis of the experimental data.  It is also a PDF in the current directory:  R. Sayag and M. G. Worster, 2013. _Axisymmetric gravity currents of power-law fluids over a rigid horizontal surface_, J. Fluid Mech. 716, doi:10.1017/jfm.2012.545.  They do not use the Halfar solution, but they do give a very clear derivation of a similarity solution very close to it.  Their solution is one which closely-matches an achievable laboratory experiment.

Read up!


needed tools
------------

* Matlab or Octave.

* Any additional tools you need to post-process results and produce a short talk at the end of the school.


questions
---------

* FIXME

goals/mileposts
---------------

* Fit the data for each Q to the rN(t) curve from the similarity solution.  Do you get the same n value(s)?

* Generate the first part of your presentation early.  I.e. slides stating the model, the parameters, the goals, and the analysis.  (You have to do this anyway, and doing it early helps reveal good questions early; that is the goal!)


recipes
-------

View the data we have from Sayag and Worster:

        $ octave    # or matlab
        >> plotSWdata

