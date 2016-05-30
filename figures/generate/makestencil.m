function makestencil
% MAKESTENCIL   Plot the kind of stencil used for illustrating finite
% difference schemes.  Works under Matlab better than Octave, especially
% because saving the figures sucks under Octave.

msize=40;

basic(1) % stencil for non-constant diffusivity
hold on
plot([-0.5 0.5 0 0], [0 0 -0.5 0.5], '^k', 'Markersize',msize)
plot([-1 0 1 0 0], [0 0 0 -1 1], 'dk', 'Markersize',msize)
hold off
% to create .eps in Octave:
%print -depsc2 ../pdffigs/diffstencil.eps

basic(2) % stencil for Mahaffy
hold on
plot([0.5], [0], '^k', 'Markersize',msize)
plot([0 1], [0 0], 'sk', 'Markersize',msize)
plot([0 1 0 1], [-1 -1 1 1], 'dk', 'Markersize',msize)
hold off
%print -depsc2 ../pdffigs/mahaffystencil.eps

basic(3) % stencil for explicit 2D heat equation
hold on
plot([-1 0 1 0 0], [0 0 0 -1 1], 'dk', 'Markersize',msize)
hold off
%print -depsc2 ../pdffigs/exp2dstencil.eps

   function basic(i)
   figure(i)
   lwidth = 8.0;

   % make the grid
   x=-1.3:0.1:1.3;
   plot(x,[-ones(size(x)); zeros(size(x)); ones(size(x))]','k',...
       'LineWidth',lwidth)
   axis([-1.7 1.7 -1.5 1.5])
   axis off
   hold on
   plot([-ones(size(x)); zeros(size(x)); ones(size(x))]',x,'k',...
       'LineWidth',lwidth)

   % label the j,k positions
   fsize = 30;
   text(-1,-1.5,'j-1','FontSize',fsize,'FontAngle','italic',...
      'HorizontalAlignment','center')
   text(0,-1.5,'j','FontSize',fsize,'FontAngle','italic',...
      'HorizontalAlignment','center')
   text(1,-1.5,'j+1','FontSize',fsize,'FontAngle','italic',...
      'HorizontalAlignment','center')
   text(-1.6,-1,'k-1','FontSize',fsize,'FontAngle','italic')
   text(-1.5,0,'k','FontSize',fsize,'FontAngle','italic')
   text(-1.62,1,'k+1','FontSize',fsize,'FontAngle','italic')
   hold off

