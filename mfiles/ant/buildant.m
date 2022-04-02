function [x,y,lat,lon,prcp,thk,topg,usrf] = buildant(doplot,filename)
% BUILDANT  Helper function for ANT, to build a Matlab/Octave model state
% for the Antarctic ice sheet, on a (default) 50km grid, from a NetCDF file.
% See comments below for how to generate a suitable NetCDF file.
% Examples:  Plots 3 fields from default file Ant50km.nc:
%   >> buildant
% This example fills the variables, but with no plot:
%   >> [x,y,lat,lon,prcp,thk,topg,usrf] = buildant(0);
% This reads from a different NetCDF file:
%   >> [x,y,lat,lon,prcp,thk,topg,usrf] = buildant(0,'Ant25km.nc');
% See also:  ANT.
%
% The preparatory steps to create Ant50km.nc require the NCO
% ("NetCDF Operators") and the download of a 104 Mb file, ALBMAP:
%   $ wget -nc http://websrv.cs.umt.edu/isis/images/4/4d/Antarctica_5km_dev1.0.nc
%   $ ncks -v lat,lon,thk,topg,usrf,acca -d x1,,,10 -d y1,,,10 \
%        Antarctica_5km_dev1.0.nc Ant50km.nc
% Higher resolution sampling:
%   $ ncks -v lat,lon,thk,topg,usrf,acca -d x1,,,5 -d y1,,,5 \
%        Antarctica_5km_dev1.0.nc Ant25km.nc
%   $ ncks -v lat,lon,thk,topg,usrf,acca -d x1,,,2 -d y1,,,2 \
%        Antarctica_5km_dev1.0.nc Ant10km.nc

if nargin < 1, doplot = 1; end
if nargin < 2, filename = 'Ant50km.nc'; end

disp(['reading variables x,y,lat,lon,acca,thk,topg,usrf from NetCDF file ' filename])

S = netcdf(filename);  % reads NetCDF file into a large structure

% the order of variables seems unclear and special to the file!
prcp = squeeze(double(S.VarArray(1).Data));
lat = squeeze(double(S.VarArray(2).Data));
lon = squeeze(double(S.VarArray(3).Data));
thk = squeeze(double(S.VarArray(5).Data));
topg = squeeze(double(S.VarArray(7).Data));
usrf = squeeze(double(S.VarArray(8).Data));
x = double(S.VarArray(9).Data);
y = double(S.VarArray(10).Data);

if doplot==0, return; end

disp('plotting 3 fields')
figure(1)
surf(x/1000,y/1000,usrf), shading('flat'), view(2), axis square
xlabel('x (km)'), ylabel('y (km)'), title('surface elevation "usrf"  (m)'), colorbar
figure(2)
surf(x/1000,y/1000,topg), shading('flat'), view(2), axis square
xlabel('x (km)'), ylabel('y (km)'), title('bed elevation "topg"  (m)'), colorbar
figure(3)
surf(x/1000,y/1000,prcp), shading('flat'), view(2), axis square
xlabel('x (km)'), ylabel('y (km)'), title('precipitation "prcp"  (m a-1)'), colorbar
