function [t vol] = ant(filename,doplot,E)
% ANT  Simulate Antarctic ice sheet flow using BUILDANT to
% extract data from re-gridded SeaRISE-Antarctic data.  See
% Ant50km.nc for metadata.
% Example: Using Ant50km.nc data for 40 ka run:
%   >> addpath('../')   % so Matlab/Octave can find other codes
%   >> ant
% Using different gridded data and different enhancement factor
% E=5.0, and saving volume time series for reload:
%   >> [t vol] = ant('Ant25km.nc',1,5.0);
%   >> save -v7 vol25km.mat t vol
% (Retrieve the data with "load vol25km.mat".)
% Run times of 9/8/2012 runs with Octave on bueler-lemur:
%   50 km = 373 sec,  25 km = 1829 sec,  20 km = 12065 sec.
% Calls:  BUILDANT, SIAGENERAL, DIFFUSION

if nargin < 1, filename = 'Ant50km.nc'; end
if nargin < 2, doplot = 1; end
if nargin < 3, E = 3; end  % default enhancement factor

% read input data from NetCDF; no plot
[x,y,lat,lon,prcp,thk,topg,usrf] = buildant(0,filename);

% grid info
Lx = (max(x) - min(x)) / 2;    Ly = (max(y) - min(y)) / 2;
J = length(x) - 1;    K = length(y) - 1;
dx = 2 * Lx / J;    dy = 2 * Ly / K;

fprintf('summary of input data:\n')
fprintf('  J = %d,  K = %d\n',J,K)
fprintf('  dx = %.3f km,  dy = %.3f km\n',dx/1000.0,dy/1000.0)
fprintf('  thickness     [min,max] = [%8.2f,%8.2f] m\n',...
  min(min(thk)),max(max(thk)))
fprintf('  bed elevation [min,max] = [%8.2f,%8.2f] m\n',...
  min(min(topg)),max(max(topg)))
fprintf('  precipitation [min,max] = [%8.5f,%8.5f] m a-1\n',...
  min(min(prcp)),max(max(prcp)))

% run-time and time-step (in years)
secpera = 31556926;
deltata = 1.0;
tblocka = 500.0;
NN = 80;  % number of blocks of length tblock
tfa = NN * tblocka;
t = (0:NN) * tblocka * secpera;
fprintf('doing run of %.3f ka total, in blocks of %.3f a,\n',...
        tfa/1000.0,tblocka)
fprintf('  with max time-steps of %.3f a ...\n  running:\n',deltata)

% fix units and parameters for actual run
deltat = deltata * secpera;  % convert to seconds
M = prcp / secpera;  % alternatively: M = zeros(size(thk));
if any(any(thk<0)), error('negative input thicknesses detected'), end
A = E * 1.0e-16 / secpera;

% solve SIA by doing blocks of tf time, and reporting volume
H = thk;
vol = printvolume(0.0,dx,dy,H);
rho = 910.0;    rhow = 1028.0;
hinit = getsurface(H,topg,rho,rhow);
for k = 1:NN
  [H,hfinal,dtlist] = siageneral(Lx,Ly,J,K,H,deltat,tblocka*secpera,topg,M,A);
  if any(any(H<0)), error('negative output thicknesses detected'), end
  vol = [vol printvolume(k*tblocka,dx,dy,H)];
end

if doplot==0, return; end

figure(1)
imagesc(x/1000,y/1000,flipud(hinit),[0, 4000]), axis square, colorbar
xlabel('x  (km)','fontsize',14), ylabel('y  (km)','fontsize',14)
%title('initial surface elevation')
%print -dpng antinitial.png

figure(2)
imagesc(x/1000,y/1000,flipud(hfinal),[0, 4000]), axis square, colorbar
xlabel('x  (km)','fontsize',14), ylabel('y  (km)','fontsize',14)
%title('final surface elevation')
%print -dpng antfinal.png

figure(3)
imagesc(x/1000,y/1000,flipud(H-thk),[-1000, 1000]), axis square, colorbar
xlabel x, ylabel y, title('thickness change')
%print -dpng antthickchange.png

figure(4)
plot(t/secpera,vol/(1.0e6*1.0e9),'o-','markersize',11,'linewidth',2)
xlabel('t  (a)','fontsize',14)
ylabel('volume  (10^6 km^3)','fontsize',14)
grid on
%print -dpdf antvol.pdf

  function vol = printvolume(timea,dx,dy,thk)
    vol = sum(sum(thk)) * dx * dy;
    fprintf('  ice volume at time %7.3fka  is %.4e km^3\n',timea/1000.0,vol/1.0e9)
  end

  function h = getsurface(H,b,rho,rhow)
    % GETSURFACE  Helper function to build surface elevations carefully,
    % according to grounded or floating.
    f = rho / rhow;                     % fraction of floating ice below surface
    h = H + b;                          % only valid where H > 0 and grounded
    h(H <= 0) = max(b(H <= 0),0.0);     % if no ice
    floating = (H > 0) & (b < - f * H); % points where flotation criterion
    h(floating) = (1-f) * H(floating);  %   ... is applied
  end
end
