function [u,u0] = ssaflowline(p,J,H,b,ug,initchoice)
% SSAFLOWLINE  Computes the velocity from the SSA in a flow-line case, with
%   left-end Dirichlet condition  u=ug  and right-end Neumann condition from
%   calving-front condition.
% form:
%   [u,u0] = ssaflowline(p,J,H,b,uleft,initchoice)
% all 6 input arguments are required:
%   p  = struct with parameters A,C,g,L,m,n,rho,rhow,secpera
%   J  = number of grid subintervals
%   H  = thickness (m), a length J+1 column vector
%   b  = bed elevation (m), same size as H
%   ug = velocity boundary value at x=0 end of domain
%   initchoice = 1,2 chooses initial guess scheme;
%                use 1 for ice shelves and 2 for ice streams
% outputs:
%   u = the velocity solution (m s-1), a length J+1 column vector
%   u0 = the velocity initial guess (m s-1), same size as u
% Does "Picard" iteration to solve nonlinear SSA problem.
% Calls:  FLOWLINE to solve inner linear PDE boundary value problem,
%         SSAINIT
% Example:  see TESTSHELF

if nargin ~= 6, error('exactly 6 input arguments required'), end

dx = p.L / J;  x = (0:dx:p.L)';
xstag = (dx/2:dx:p.L+dx/2)';  % yes, it goes past end

% create parts of PDE problem to solve; see flowline.m
% coefficient for dragging:
alpha = p.C * ones(size(x));  % applies to linear till only; change for nonlinear
h = H + b;
hx = regslope(dx,h);
beta = p.rho * p.g * H .* hx;  % driving stress becomes right-hand side
% value in calving front stress boundary condition:
gamma = ( 0.25 * p.A^(1/p.n) * (1 - p.rho/p.rhow) *...
          p.rho * p.g * H(end) )^p.n;

u0 = ssainit(p,x,beta,gamma,initchoice);  u = u0;

% "outer" iteration for solution of nonlinear equation;
%   "Picard" iteration
Hstag = stagav(H);
tol = 1.0e-14;  % m s-1; = 0.000316 mm a-1 velocity error
eps_reg = (1.0 / p.secpera) / p.L;  % strain rate of 1 m/a over length of shelf
maxdiff = Inf;  W = zeros(J+1,1);  iter = 0;
while maxdiff > tol
  % find coefficient W(x) on staggered grid using "old" u
  uxstag = stagslope(dx,u);
  sqr_ux_reg = uxstag.^2 + eps_reg^2; % regularize to avoid division by zero
  W(1:J) = 2 * p.A^(-1/p.n) * Hstag .* sqr_ux_reg.^(((1/p.n)-1)/2.0);
  W(J+1) = W(J);  % past end:  local constancy of effective viscosity
                  %            times thickness

  % solve problem for this W
  unew = flowline(p.L,J,gamma,W,alpha,beta,ug);
  maxdiff = max(abs(unew-u));
  u = unew;
  iter = iter + 1;
  fprintf('.')
end
fprintf('\nSSA solver did %d Picard iterations on dx = %.3f km grid\n',iter,dx/1000.0)

%STRIPFROMHERE <- only relevant to what is shown in slides

  function fav = stagav(f)
  % average regular grid values onto staggered grid
  fav = 0.5 * (f(1:end-1) + f(2:end));

  function slope = regslope(dx,f)
  % compute regular grid values of slope  f_x = f'
  % if length(h)=J+1, returns vector with same length J+1
  J = length(f) - 1;
  slope = [(f(2)-f(1))/dx; (f(3:J+1)-f(1:J-1))/(2*dx); (f(J+1)-f(J))/dx];

  function slope = stagslope(dx,f)
  % compute staggered grid values of slope  f_x = f'
  % if length(h)=J+1, returns vector with length J
  slope = (f(2:end) - f(1:end-1)) / dx;
