function [averr,maxerr] = testshelf(J,L)
% TESTSHELF  calls ssaflowline.m to compute a numerical solution, calls 
%   exactshelf.m to compute the corresponding exact solution, and reports the
%   average and maximum velocity errors
% form:
%   maxerr = testshelf(J,L)
% where:
%   J  = number of grid points
%   L  = length of ice shelf (m)
% Example: 
%   >> testshelf(25,200e3);        % show figure, coarse grid
%   >> err = testshelf(400,200e3)  % no figure, return vel error in m/s
% Called by:  SHELFCONV
% Calls:  SSAFLOWLINE

% use a structure for physical parameters (values are from MISMIP)
param = struct('secpera',31556926,...
               'n',3.0,...       % Glen flow law exponent
               'rho',900.0,...   % kg m^-3
               'rhow',1000.0,... % kg m^-3
               'g',9.8,...       % m s-2
               'A',1.4579e-25);  % s^-1 Pa^-3; used by MacAyeal et al (1996) for
                                 % EISMINT-Ross; corresponds to B=1.9e8 Pa s^1/3

% default shelf has length 200 km, accumulation of 30 cm/a and
%   thickness 500 m and velocity 50 m/a at grounding line
Hg = 500;  ug = 50 / param.secpera;
M0 = 0.3 / param.secpera;
if nargin<2, param.L = 200e3; else, param.L = L; end
if nargin<1, J = 300; end

dx = param.L / J;   x = (0:dx:param.L)';

% get exact solution
[uexact,H] = exactshelf(x,param.L,M0,Hg,ug);

% get numerical
r = param.rho / param.rhow;  b = - r * H;  % flotation criterion
param.C = 0;  % no basal drag for a shelf
[unum,u0] = ssaflowline(param,J,H,b,ug,1);
averr = sum(abs(unum-uexact)) / (J+1);
maxerr = max(abs(unum-uexact));

if nargout==0  % plot, if not asked for error
  plot(x/1000.0,u0 * param.secpera,'b-','LineWidth',3.0),  hold on
  plot(x/1000.0,uexact * param.secpera,'r-','LineWidth',3.0)
  plot(x/1000.0,unum * param.secpera,'go','markersize',12,'LineWidth',2.0),  hold off
  xlabel('x  (km)')%,  ylabel('velocity  (m/a)')
  %title(sprintf('case with  dx = %.0f km  has  (max error) = %.3f  m a^{-1}',...
  %               dx/1000.0,maxerr * param.secpera))
end

