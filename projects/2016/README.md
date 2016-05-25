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
        ln -s ~/repos/sia-fve/petsc/figsmahaffy.py
        ln -s ~/petsc/bin/PetscBinaryIO.py
        ln -s ~/petsc/bin/petsc_conf.py

Generate and view example bed:

        ./bedrand.py --plotbed

Generate and save example bed:

        ./bedrand.py --plotbed -o low.dat

Run a low-resolution example using the bed:

        mkdir testlow
        ./mahaffy -mah_read low.dat -mah_showdata -draw_pause 2 -mah_cmbmodel -cmb_ela 2500.0 -cmb_lapse 0.002 -cs_D0 1.0 -mah_dump testlow/

Look at result:

        cd testlow/
        ../figsmahaffy.py --observed

Run a higher-resolution example, in parallel:

        ./bedrand.py --plotbed -Nx 100 -Ny 100 -o high.dat
        mkdir testhigh
        mpiexec -n 4 ./mahaffy -mah_read high.dat -mah_showdata -draw_pause 2 -mah_cmbmodel -cmb_ela 2500.0 -cmb_lapse 0.002 -cs_D0 1.0 -mah_dtrecovery 10 -snes_max_it 400 -pc_type asm -sub_pc_type lu -mah_dump testhigh/
        cd testhigh/
        ../figsmahaffy.py --observed


Ice rheology from symmetric flows
---------------------------------

FIXME



