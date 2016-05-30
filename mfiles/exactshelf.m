function [u,H] = exactshelf(x,L,M0,Hg,ug)
% EXACTSHELF  computes the exact thickness and velocity in a steady ice shelf,
%   with grounding line at xg = 0, length L, constant surface mass balance M0,
%   and given thickness (Hg) and velocity (ug) at grounding line
% form:
%   [u,H] = exactshelf(x,L,M0,Hg,ug)
% where:
%   x  = grid of x-values
%   L  = length of ice shelf (m)
%   M0 = surface mass balance rate (m/s)
%   Hg = thickness at grounding line (m)
%   ug = velocity at grounding line (m/s)
% example:
%   >> L=200e3;  x=linspace(0,L,1000);                   
%   >> [u,H] = exactshelf(x,L);   % allow defaults for M0,Hg,ug
%   >> h = 0.1 * H; b = -0.9 * H;  % flotation surface elevation and base elev
%   >> figure(1), plot(x/1000,h,x/1000,b), xlabel('x  (km)'), ylabel('elev (m)')
%   >> figure(2), plot(x/1000,u*3.1556926e7), xlabel('x  (km)'), ylabel('u (m/a)')


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
if nargin<5, ug = 50 / param.secpera; end
if nargin<4, Hg = 500; end
if nargin<3, M0 = 0.3 / param.secpera; end
if nargin<2, L = 200e3; end
if nargin<1, error('x values required as first argument'), end

% see lecture notes for analysis leading to exact ice shelf shape
n = param.n;
r = param.rho/param.rhow;
Cs = param.A * (0.25 * param.rho * param.g * (1-r))^n;
qg = ug * Hg;

u = zeros(size(x)); H = u;
in = (x<=L) & (x>=0);
u(in) = ( ug^(n+1) + (Cs/M0) * ((M0 * x(in) + qg).^(n+1) - qg^(n+1)) ).^(1/(n+1));
H(in) = (M0 * x(in) + qg) ./ u(in);

