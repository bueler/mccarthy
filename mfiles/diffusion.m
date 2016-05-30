function [T,dtav] = diffusion(Lx,Ly,J,K,Dup,Ddown,Dright,Dleft,T0,tf,F,b)
% DIFFUSION  Adaptive explicit method for diffusion equation
%   T_t = F + div (D grad (T + b))
% on rectangle (-Lx,Lx) x (-Ly,Ly) with initial condition T0 and
% functions F(x,y), D(x,y), b(x,y).  The boundary condition is that
% the T values at the edge of the domain (i.e. T(1,.), T(J+1,.), T(.,1),
% T(.,K+1)) are all held fixed at their initial values supplied in T0.
% That is, the boundary condition is a Dirichlet condition from the
% initial values.
% Usage:
%   T = diffusion(Lx,Ly,J,K,Dup,Ddown,Dright,Dleft,T0,tf,F)
% where
%   T     = approximate solution at tf
%   Lx,Ly = half-widths of rectangular domain
%   J,K   = number of subintervals in x,y directions, resp.
%   D*    = (J-1) x (K-1) arrays with diffusivities for "staggered" grid
%   T0    = (J+1) x (K+1) array with initial values on regular grid
%   tf    = final time
%   F     <-- OPTIONAL ARGUMENT
%         = (J+1) x (K+1) array with heat source on regular grid
%   b     <-- OPTIONAL ARGUMENT
%         = (J+1) x (K+1) array with offset before gradient
% Note: There is no error checking on sizes of D*, T0, F, b arrays.
% Note: There is no text output.  Restore commented 'fprintf' if desired.
% Note: The input diffusivities could be time-dependent, but they are
%   time-independent in this simplified implementation.  Call-back could
%   implement time-dependent diffusivity D(t,x,y).
% Example: Compare this result to HEATADAPT:
%   >> J = 50;  K = 50;  D = ones(J-1,K-1);
%   >> [x,y] = ndgrid(-1:2/J:1,-1:2/K:1);
%   >> T0 = exp(-30*(x.*x + y.*y));
%   >> T = diffusion(1.0,1.0,J,K,D,D,D,D,T0,0.05);
%   >> surf(x,y,T), shading('interp'), xlabel x, ylabel y
% Called by: SIAFLAT, SIAGENERAL.

% spatial grid and initial condition
dx = 2 * Lx / J;    dy = 2 * Ly / K;
[x,y] = ndgrid(-Lx:dx:Lx, -Ly:dy:Ly); % (J+1) x (K+1) grid in x,y plane
T = T0;
if nargin < 11, F = zeros(size(T0)); end  % for heat source term
if nargin < 12, b = zeros(size(T0)); end  % allows use for nonflat-bed SIA case

%fprintf('  doing explicit steps adaptively on 0.0 < t < %.3f\n',tf)
t = 0.0;    count = 0;
while t < tf
   % stability condition gives time-step restriction
   maxD = [max(max(Dup)) max(max(Ddown)) max(max(Dleft)) max(max(Dright))];
   maxD = max(maxD);  % scalar maximum of D
   if maxD <= 0.0  % e.g. happens with zero thickness ice sheets
     dt = tf - t;
   else
     dt0 = 0.25 * min(dx,dy)^2 / maxD;
     dt = min(dt0, tf - t);  % do not go past tf
   end
   mu_x = dt / (dx*dx);    mu_y = dt / (dy*dy);
   Tb = T + b;
   T(2:J,2:K) = T(2:J,2:K) + ...
       mu_y * Dup    .* ( Tb(2:J,3:K+1) - Tb(2:J,2:K)   ) - ...
       mu_y * Ddown  .* ( Tb(2:J,2:K)   - Tb(2:J,1:K-1) ) + ...
       mu_x * Dright .* ( Tb(3:J+1,2:K) - Tb(2:J,2:K)   ) - ...
       mu_x * Dleft  .* ( Tb(2:J,2:K)   - Tb(1:J-1,2:K) );
   T = T + F * dt;
   t = t + dt;    count = count + 1;
   %fprintf('.')
end
dtav = tf / count;
%fprintf('\n  completed N = %d steps, average dt = %.7f\n',count,tf/count)
%surf(x,y,T),  shading('interp'),  xlabel x,  ylabel y
