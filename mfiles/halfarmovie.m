% HALFARMOVIE   Show Halfar (1983) solution by frames.  Calls HALFAR.
% (Uncomment print line to create PNG figures which form movie shown in the
% lecture.  Naming of files allows "animate" LaTeX package to read from anim/.)

L = 1200e3;  dx = L/100;  [x,y] = meshgrid(-L:dx:L,-L:dx:L);
count = 0;
for tyears = 10.0.^(-0.5:0.25:6)
  H = halfar(tyears * 3.1556926e7,x,y);
  figure(1), clf
  surf(x/1000,y/1000,H),  shading('interp')
  xlabel('x  (km)'),  ylabel('y  (km)'),  zlabel('H  (m)')
  axis tight, view(3)
  axis([-L/1000 L/1000 -L/1000 L/1000 0 8000])
  %print(sprintf('../anim/halfar%d.png',count), '-dpng', '-r100')
  count = count + 1;
  sleep(0.5)
end

