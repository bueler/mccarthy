function [H,h,dtlist] = siageneral(Lx,Ly,J,K,H0,deltat,tf,b,M,A)
% SIAGENERAL  Numerical solution of isothermal n=3 SIA
%   H_t = M + div (D grad (H + b))
% where  H = (ice thickness), the ice surface elevation is h = H + b,
% M = (surface mass balance), and D is the diffusivity.  We allow
% with general bed elevation b(x,y), general surface mass balance M(x,y),
% and nonnegative thickness rule.  Any floating ice is immediately calved.
% Compare SIAFLAT.
% Usage:
%   [H,h,dtlist] = siageneral(Lx,Ly,J,K,H0,deltat,tf,b,M,A)
% where
%   H      = numerical approx of thickness at final time
%   h      = numerical approx of surface elevation at final time
%   dtlist = list of time steps used adaptively in diffusion.m
%   Lx,Ly  = half lengths of rectangle in x,y directions
%   J,K    = number of subintervals in x,y directions
%   H0     = initial thickness, a (J+1)x(K+1) array
%   deltat = major time step
%   tf     = final time
%   b      = bed elevation, a (J+1)x(K+1) array
%   M      = surface mass balance, a (J+1)x(K+1) array
%   A      = ice softness
% Calls: DIFFUSION, which does adaptive explicit time-stepping within the
% major time step deltat.
% Example: See ANT.

% constants
g = 9.81;    rho = 910.0;    rhow = 1028.0;
Gamma  = 2 * A * (rho * g)^3 / 5; % see Bueler et al (2005)
f = rho / rhow;                   % fraction of floating ice below surface

dx = 2 * Lx / J;   dy = 2 * Ly / K;
N = ceil(tf / deltat);   deltat = tf / N;
j  = 2:J;   k = 2:K;  % interior indices
nk = 3:K+1;   sk = 1:K-1;   ej = 3:J+1;   wj = 1:J-1; % north,south,east,west

t = 0;   dtlist = [];
H = H0;
for n=1:N
  % staggered grid thicknesses
  Hup  = 0.5 * ( H(j,nk) + H(j,k) ); % up
  Hdn  = 0.5 * ( H(j,k) + H(j,sk) ); % down
  Hrt  = 0.5 * ( H(ej,k) + H(j,k) ); % right
  Hlt  = 0.5 * ( H(j,k) + H(wj,k) ); % left
  h = max(0.0,H + b); % update surface; uses sea level of z=0
  % staggered grid value of |grad h|^2 = "alpha^2"
  a2up = (h(ej,nk) + h(ej,k) - h(wj,nk) - h(wj,k)).^2 / (4*dx)^2 + ...
         (h(j,nk) - h(j,k)).^2 / dy^2;
  a2dn = (h(ej,k) + h(ej,sk) - h(wj,k) - h(wj,sk)).^2 / (4*dx)^2 + ...
         (h(j,k) - h(j,sk)).^2 / dy^2;
  a2rt = (h(ej,k) - h(j,k)).^2 / dx^2 + ...
         (h(ej,nk) + h(j,nk) - h(ej,sk) - h(j,sk)).^2 / (4*dy)^2;
  a2lt = (h(j,k) - h(wj,k)).^2 / dx^2 + ...
         (h(wj,nk) + h(j,nk) - h(wj,sk) - h(j,sk)).^2 / (4*dy)^2;     
  % Mahaffy's staggered grid diffusivity: D = Gamma H^{n+2} |grad h|^{n-1}
  Dup  = Gamma * Hup.^5 .* a2up;
  Ddn  = Gamma * Hdn.^5 .* a2dn;
  Drt  = Gamma * Hrt.^5 .* a2rt;
  Dlt  = Gamma * Hlt.^5 .* a2lt;
  % call *adaptive* diffusion() to time step H
  [H,dtadapt] = diffusion(Lx,Ly,J,K,Dup,Ddn,Drt,Dlt,H,deltat,M,b);
  H = max(H,0.0);       % enforce nonnegative thickness (H->0 can occur where M<0)
  calvehere = (b < - f * H);  % floating
  H(calvehere) = 0.0;   % calving occurs anywhere the ice is floating
  t = t + deltat;
  dtlist = [dtlist dtadapt];
end
h = max(0.0,H + b); % finalize surface
