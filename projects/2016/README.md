projects/2016/
==============

These are notes and data and codes for my suggested projects for the McCarthy
summer school in 2016.

Glaciating random terrain
-------------------------

Build PETSc code

        (cd ~/repos/sia-fve/petsc/ && make mahaffy)

Generate links:

        ln -s ~/repos/sia-fve/petsc/mahaffy
        ln -s ~/petsc/bin/PetscBinaryIO.py
        ln -s ~/petsc/bin/petsc_conf.py

Generate and view example bed:

        ./bedrand.py --plotbed

Generate and save example bed:

        ./bedrand.py --plotbed -o foo.dat

FIXME converges trivially:

        ./mahaffy -mah_read foo.dat -mah_showdata -draw_pause 2

Run an example  FIXME does not work yet:

        ./mahaffy -mah_read foo.dat -mah_showdata -draw_pause 2 -mah_cmbmodel -cmb_ela 3000.0 -cmb_lapse 0.001




Ice rheology from symmetric flows
---------------------------------

FIXME



