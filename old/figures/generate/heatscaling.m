% plots scaled versions of fundamental solution to heat equation

x = -6:.01:6;
y = exp(-x.^2/4);
JJ = 1:4;

tt = 10.^((JJ-2)/3)
tmalpha = tt.^(-1);
tbeta = tt.^(1/2);

figure(1), clf
for j=JJ
  subplot(1,length(JJ),j)
  xs = tbeta(j) * x;       % x scaled
  ys = tmalpha(j) * y;     % y scaled
  % plot graph and tight box around it
  plot(xs,ys,'k','LineWidth',3.0)
  hold on
  mid = floor((length(x)+1)/2);
  plot([xs(1) xs(end) xs(end) xs(1) xs(1)],[ys(1) ys(1) ys(mid) ys(mid) ys(1)])
  hold off
  axis([-15 15 0 2.5])
  axis off
end

print -dpdf heatscaling.pdf

