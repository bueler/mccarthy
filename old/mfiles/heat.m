function T = heat(D,J,K,dt,N)
% HEAT  The explicit method for the 2D heat equation
%   T_t = D (T_xx + T_yy)
% on -1 < x < 1, -1 < y < 1, for 0 < t < N * dt
% Uses fixed time steps; compare HEATADAPT.  Uses a specific
% gaussian initial condition.
% Usage:
%   T = heat(D,J,K,dt,N)
% where
%   T   = approximate solution at tf = N * dt
%   D   = diffusivity coeff
%   J,K = number of subintervals in x,y directions, resp.
%   dt  = fixed time step
%   N   = number of time steps
% Examples:
%   >> heat(1.0,30,30,0.001,0);    % just show initial condition
%   >> heat(1.0,30,30,0.001,20);   % final time 0.02
%   >> heat(1.0,30,30,0.004,5);    % final time 0.02; huh?

% setup spatial grid and initial condition:
dx = 2 / J;    dy = 2 / K;
[x,y] = meshgrid(-1:dx:1, -1:dy:1); % (J+1) x (K+1) grid in x,y plane
T = exp(-30*(x.*x + y.*y));
fprintf('  doing N = %d steps of dt = %.5f for 0.0 < t < %.3f\n',N,dt,N*dt)
fprintf('  nu = 1 * dt / dx^2 = %.5f\n',dt/dx^2)

% do explicit time steps
mu_x = dt * D / (dx*dx);
mu_y = dt * D / (dy*dy);
for n=1:N
   % colon notation substitutes for double loop over interior grid points
   T(2:J,2:K) = T(2:J,2:K) + ...
       mu_x * ( T(3:J+1,2:K) - 2 * T(2:J,2:K) + T(1:J-1,2:K) ) + ...
       mu_y * ( T(2:J,3:K+1) - 2 * T(2:J,2:K) + T(2:J,1:K-1) );
end

% show solution
surf(x,y,T),  shading('interp'),  xlabel x,  ylabel y
