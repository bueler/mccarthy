# mfiles/ant/

A 50 km resolution NetCDF data set is already prepared.  This is `Ant50km.nc`.  To run the default 40 ka Antarctica simulation using this 50 km grid do
To run the codes in the `ant/` subdirectory, do this:

        >> addpath('../')   % so Matlab/Octave can find the other codes
        >> ant;             % run ant.m

The preparatory steps to create Ant50km.nc require the [NCO ("NetCDF Operators")](http://nco.sourceforge.net/) and the download of a 104 Mb file, ALBMAP v1:
        $ wget -nc http://websrv.cs.umt.edu/isis/images/4/4d/Antarctica_5km_dev1.0.nc
See `albmap.bib` for the citation to this data.

Next one needs to sub-sample using the `ncks` operator:
        $ ncks -v lat,lon,thk,topg,usrf,acca -d x1,,,10 -d y1,,,10  Antarctica_5km_dev1.0.nc Ant50km.nc

For higher-resolution sub-sampling:
        $ ncks -v lat,lon,thk,topg,usrf,acca -d x1,,,5 -d y1,,,5  Antarctica_5km_dev1.0.nc Ant25km.nc
        $ ncks -v lat,lon,thk,topg,usrf,acca -d x1,,,2 -d y1,,,2  Antarctica_5km_dev1.0.nc Ant10km.nc

Now run `ant.m` with the filename and additional options:
        >> [t vol] = ant('Ant25km.nc',1,5.0);
Do
        >> help ant
for more help.