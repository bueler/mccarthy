% HALFARMOVIE   Show Halfar (1983) solution by frames.  Calls HALFAR.
% (Uncomment print line to generate PNG figures for animation in slides.)

L = 1200e3;  dx = L/100;  [x,y] = meshgrid(-L:dx:L,-L:dx:L);
count = 0;
for tyears = [0.3 1 3 10 30 100 300 1000 3000 10000 30000 100000 300000 1000000]
  H = halfar(tyears * 3.1556926e7,x,y);
  figure(1), clf
  surf(x/1000,y/1000,H),  shading('interp')
  xlabel('x  (km)'),  ylabel('y  (km)'),  zlabel('H  (m)')
  axis tight, view(3)
  axis([-L/1000 L/1000 -L/1000 L/1000 0 8000])
  if tyears < 1
      timelabel = sprintf('t=%.2f years',tyears);
  else
      timelabel = sprintf('t=%d years',tyears);
  end
  text(-900,-500,9000,timelabel,'fontsize',16)
  %print(sprintf('../slides/anim/halfar-%d.png',count), '-dpng', '-r100')
  count = count + 1;
  pause(0.5)
end

