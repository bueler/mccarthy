function maxerr = testflowline(test,J)
% TESTFLOWLINE  tests flowline.m by setting up and solving problems:
%   test=1:  trivial straight line solution
%   test=2:  manufactured solution with W(x) constant
%   test=3:  manufactured solution with W(x) non-constant
% form:
%   maxerr = testflowline(test,J)
% where:
%   test = explained above
%   J    = grid has J subintervals
% example: see convanalysis.m

if nargin<2, J=300; end
if nargin<1, test=2; end

L = 1;
dx = L / J;
x = (0:dx:L)';
xstag = (dx/2:dx:L+dx/2)';  % last staggered point is past end
if test == 1
  W = ones(size(xstag));
  alpha = zeros(size(x));
  beta = zeros(size(x));
  gamma = 1 / L;
  uexact = x / L;
elseif test == 2
  W = 2 * ones(size(xstag));
  alpha = ones(size(x));
  uexact = sin(x);
  gamma = cos(1);  % exact derivative at end
  % manufacture beta by   beta := (W u_x)_x - alpha u
  beta = - 3 * sin(x);
elseif test == 3
  W = 2 - xstag;
  alpha = ones(size(x));
  uexact = sin(x);
  gamma = cos(1);  % exact derivative at end
  % manufacture beta by   beta := (W u_x)_x - alpha u
  beta = - cos(x) - (3 - x) .* sin(x);
else, error('only test==1,2,3 allowed'), end

% run flowline
unum = flowline(L,J,gamma,W,alpha,beta,0);
maxerr = max(abs(unum-uexact));

if nargout==0
  plot(x,uexact,x,unum,'o');
  fprintf('maximum absolute error = %.4e\n',maxerr);
end  

