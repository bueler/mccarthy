function heatstencils
% HEATSTENCILS   Plot the kind of stencils used for explicit, implicit,
% and Crank-Nicolson methods on 1D heat equation  u_t = D u_xx.

msize=20;

basic(1) % stencil for explicit
hold on
plot([-1 0 1], [0 0 0], 'dk', 'Markersize',msize)
plot([0], [1], 'sk', 'Markersize',msize)
hold off
% to create .eps in Octave:
%print -depsc2 ../pdffigs/expstencil.eps

return

basic(2) % stencil for implicit
hold on
plot([0], [0], 'dk', 'Markersize',msize)
plot([-1 0 1], [1 1 1], 'sk', 'Markersize',msize)
hold off
%print -depsc2 ../pdffigs/impstencil.eps

basic(3) % stencil for Crank-Nicolson
hold on
plot([-1 0 1], [0 0 0], 'dk', 'Markersize',msize)
plot([-1 0 1], [1 1 1], 'sk', 'Markersize',msize)
hold off
%print -depsc2 ../pdffigs/cnstencil.eps

   function basic(i)

   lwidth = 8.0;
   fsize = 20;

   figure(i), clf

   % make the grid
   x=-1.3:0.1:1.3;
   plot(x,[zeros(size(x)); ones(size(x))]','k',...
       'LineWidth',lwidth)
   axis([-1.7 1.7 -0.5 1.5])
   axis off
   hold on
   y=-0.3:0.1:1.3;
   plot([-ones(size(y)); zeros(size(y)); ones(size(y))]',y,'k',...
       'LineWidth',lwidth)

   % label the j,k positions
   text(-1,-0.4,'j-1','FontSize',fsize,'FontAngle','italic',...
      'HorizontalAlignment','center')
   text(0,-0.4,'j','FontSize',fsize,'FontAngle','italic',...
      'HorizontalAlignment','center')
   text(1,-0.4,'j+1','FontSize',fsize,'FontAngle','italic',...
      'HorizontalAlignment','center')
   text(-1.7,0,'n','FontSize',fsize,'FontAngle','italic')
   text(-1.7,1,'n+1','FontSize',fsize,'FontAngle','italic')
   hold off

