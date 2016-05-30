function H = halfar(t,x,y)
% HALFAR  Compute the similarity solution to the isothermal flat-bed SIA from
% Halfar (1983).  Constants H0 = 3600 m and R0 = 750 km are as in Test B in
% Bueler et al (2005).  Usage:
%   H = halfar(t,x,y)
% where
%   x,y = grid; can be matrices, of same size; in meters
%   t   = time in seconds (scalar)
% Example:
%   >> L = 1000e3;  dx = L/100;  [x,y] = meshgrid(-L:dx:L,-L:dx:L);
%   >> H = halfar(100.0*3.1556926e7,x,y);
%   >> surf(x,y,H),  shading('interp'),  xlabel x,  ylabel y

g = 9.81;     % constants in SI units
rho = 910.0;
secpera = 31556926;
n = 3;
A = 1.0e-16/secpera; 
Gamma  = 2 * A * (rho * g)^3 / 5; % see Bueler et al (2005)

H0 = 3600;
R0 = 750e3;
alpha = 1/9;
beta = 1/18;
% for t0, see equation (9) in Bueler et al (2005); result is 422.45 a:
t0 = (beta/Gamma) * (7/4)^3 * (R0^4/H0^7);

r = sqrt(x.*x + y.*y);
r = r / R0;
t = t / t0;
inside = max(0, 1 - (r / t^beta).^((n+1) / n));
H = H0 * inside.^(n / (2*n+1)) / t^alpha;
