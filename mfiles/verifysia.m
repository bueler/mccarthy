function [avthkerr,maxthkerr] = verifysia(J,dtyears)
% VERIFYSIA  Compare the Halfar (1983) similarity solution at year
% t = 20 ka to the numerical solution using SIAFLAT, for a run from
% t = 200 a to t = 20 ka.  The t = 200 a state of the Halfar solution
% is used as the initial value for the run.
% Usage:
%   verifysia(J,dtyears)
% where:
%   J = number of grid spaces in x and y directions
%   dtyears = "major" time step, during which diffusivity does not change;
%     note that diffusion() is called by siaflat.m and it does its own
%     adaptive time stepping
% Example:
%   >> verifysia;

if nargin<1, J=40; end;
if nargin<2, dtyears=10.0; end;

L = 1200e3;    dx = 2 * L / J;
[x,y] = meshgrid(-L:dx:L, -L:dx:L);

t1 = 200;    t2 = 20000;    secpera = 31556926;
H1 = halfar(t1 * secpera,x,y); % initial condition

[H2approx,dtlist] = siaflat(L,L,J,J,H1,dtyears*secpera,(t2-t1)*secpera);
% alternatively, test siageneral:
%[H2approx,h2approx,dtlist] = siageneral(L,L,J,J,H1,...
%                               dtyears*secpera,(t2-t1)*secpera,...
%                               zeros(J+1,J+1),zeros(J+1,J+1),1.0e-16/secpera);
% perhaps measure mass (volume) conservation:
%format long e
%vol1 = sum(sum(H1))*dx*dx
%vol2approx = sum(sum(H2approx))*dx*dx
%voldiff = abs(vol2approx - vol1)
%relvoldiff = voldiff / vol1

H2exact = halfar(t2 * secpera,x,y);
error = H2approx - H2exact;

fprintf('errors for %d x %d grid:\n',J,J)
avthkerr = sum(sum(abs(error)))/(J+1)^2;
fprintf('average thickness error  = %.3f\n',avthkerr)
maxthkerr = max(max(abs(error)));
fprintf('maximum thickness error  = %.3f\n',maxthkerr)
erreta = H2approx.^(8/3) - H2exact.^(8/3);
maxeta = max(max(H2exact.^(8/3)));
fprintf('rel max abs error in H^{8/3} = %.5f\n',max(max(abs(erreta / maxeta))) )

figure(1),   surf(x,y,error)
shading('flat'),   axis square,   view(2),   colorbar
title('thickness error (m)')

% side-by-side comparison of numerical and exact result:
figure(2),   subplot(121),   surf(x/1000,y/1000,H2exact),   shading('flat')
xlabel('x (km)'),   ylabel('y (km)'),   zlabel('exact thickness (m)')
subplot(122),   surf(x/1000,y/1000,H2approx),   shading('flat')
xlabel('x (km)'),   ylabel('y (km)'),   zlabel('numerical thickness (m)')

% to make figure showing adaptive time-stepping:
% figure(3), plot(dtlist / secpera,'o')
% xlabel('step'), ylabel('length of step in years')
