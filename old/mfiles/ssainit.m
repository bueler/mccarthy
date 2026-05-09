function u0 = ssainit(p,x,drivestress,gamma,initchoice)
% SSAINIT  Computes the initial guess velocity u0, used
%          in ssaflowline.m to solve the SSA flow line
%          stress balance.
% form
%   u0 = ssainit(p,x,drivestress,gamma,initchoice)
% all 5 input arguments are required:
%   p as in ssaflowline.m
%   x is (J+1) length column vector of grid locations
%   drivestress = rho g H h_x = driving stress; same size as x
%   gamma as in ssaflowline.m
%   initchoice = 1,2,3 chooses initial guess scheme;
%                use 1 for ice shelves and 2 for ice streams
% outputs:
%   u0 = the velocity initial guess (m s-1), same size as x
% Called by SSAFLOWLINE

if nargin ~= 5, error('exactly 5 input arguments required'), end

if initchoice == 1
  % linear initial guess from u_xx = 0, u(0)=0, and fraction
  %   of calving condition;  APPROPRIATE TO ICE SHELVES
  u0 = gamma * x;
elseif initchoice == 2
  % initial guess is that velocity is function of driving stress
  %   because resistance is entirely at base with no membrane
  %   resistance;  APPROPRIATE TO ICE STREAMS
  u0 = ( -(1/p.C) * drivestress ).^(1/p.m);
elseif initchoice == 3
  % linear velocity based on seen solution in one case:
  %   u = 30 + 30 * x / L 
  u0 = (30.0 / p.secpera) * ones(size(x)) +...
       (30.0 / p.secpera) * (x / p.L);
else, error('initchoice must be 1,2,3'), end

