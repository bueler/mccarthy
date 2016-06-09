
Ice rheology from symmetric flows
---------------------------------


reassurance
-----------

The point of the project is not to get you to a particular spot!  The point is to _help you learn fast from whereever you currently are_.

So please be honest with yourself about what you don't know, and please ask about that!  Another student may program better than you, or have seen more differential equations, or have more glaciers background, or whatever.  This is irrelevant to my goal, which is to help you move fast, in one intense week.

So, though I want to help you go forward in the direction of this project and topic, it must be from the things _you_ do understand into things _you_ do not understand.


description
-----------

The original description:

> PROJECT 7:  Ice rheology from symmetric flows

> STUDENTS: XXX

> ADVISOR: Ed Bueler

> DESCRIPTION:  There are simple models for the evolution of circularly-symmetric flows of ice (or other viscous fluids).  This theory, used quantitatively, allows us to convert observations of real flows, whether approximately circular (polar caps on Mars) or deliberately circular (laboratory), into estimates of the rheology of the fluid.  The project will use pencil and paper calculations, some processing of data on computer, and modest experimentation.

> SOFTWARE REQUIREMENTS: Matlab or Octave

> REQUIRED STUDENT BACKGROUND: Basic exposure to (1) differential equations and (2) Matlab or similar.


references and reading
----------------------

First, see the notes (`notes/notes.pdf`).  They cover the shallow ice approximation (SIA; equation (26) in notes) which you will want to comprehend, and the Halfar (1983) similarity solution (equation (40) in notes).

The Halfar solution is used in the following paper, a PDF of which is in the current directory: J. F. Nye, 2000. _A flow model for the polar caps of Mars_, J. Glaciol. 46 (154), 438--444.  It has a good derivation of the Halfar solution.

More of the Nye argument that the caps are not primarily CO2 ice is in this reference, for which there is also a PDF: Nye, J. F., Durham, W. B., Schenk, P. M., & Moore, J. M., 2000. _The instability of a south polar cap on Mars composed of carbon dioxide_, Icarus 144 (2), 449--455.

Also the Halfar solution is derived in the 2nd edition of van der Veen's textbook.  Namely starting on page 313 of C. J. van der Veen, 2013. _Fundamentals of Glacier Dynamics_, 2nd ed., CRC Press.  A copy of this textbook is in the lecture room.

The single most important paper to read is this one, for which I am proposing you reproduce part of the analysis of the experimental data.  It is also a PDF in the current directory:  R. Sayag and M. G. Worster, 2013. _Axisymmetric gravity currents of power-law fluids over a rigid horizontal surface_, J. Fluid Mech. 716, doi:10.1017/jfm.2012.545.  They do not use the Halfar solution, but they do give a very clear derivation of a similarity solution which is closely-related to it.  Their solution is one which closely-matches an achievable laboratory experiment.

Read up and ask lots of questions!


needed tools
------------

* Matlab or Octave.

* Any additional tools you need to post-process results and produce a short talk at the end of the school.


questions
---------

* What are the fundamental assumptions in the SIA?  To what extent do they apply to the Mars polar caps?  The constant-flux experiment in Sayag & Worster (2013)?

* How is the similarity solution in Sayag & Worster (2013) different from the Halfar solution?

* How do you fit data to a curve of the form  f(t) = C t^q  so that you can recover the best values of C and q?

* Why do Sayag & Worster (2013) NOT use all the data to determine their best fit values of n and mu?

* Suppose you did not have an exact (e.g. similarity) solution.  Could you re-do the experimental analysis without the exact solution but instead only using numerical solutions?  What are the issues?


goals/mileposts
---------------

* Be able to derive the velocity formula that applies in the SIA, following the way the notes do it by using the slab-on-slope case as though it is generic.

* Similarly, be able to derive the SIA evolution equation for the thickness, equation (26) in the notes.

* Follow and understand the derivation of the Halfar solution in Nye (2000) or van der Veen.

* Fit the data for each Q (three different values) to the rN(t) curve from the similarity solution in Sayag & Worster (2013).  Do you get the same n value(s) as they do?  Which data do they fit to?

* Generate the first part of your presentation early.  I.e. slides stating the model, the parameters, the goals, and the analysis.  In the initial draft, _make up_ the analysis; this exercise will help you understand what you want to do.  (_You have to build the presentation anyway, and doing it early helps reveal good questions early.  Good questions are the goal!_)


recipes
-------

View the data we have from Sayag and Worster, from the current directory:

        $ octave    # or matlab
        >> plotSWdata

Generate some random data, plot it, and do linear regression on it:

        >> x = 0:0.5:4
        >> y = x/2 + rand(1,9)
        >> plot(x,y,'o')
        >> p = polyfit(x,y,1)
        >> plot(x,y,'o',x,p(1)*x+p(2),'--')

