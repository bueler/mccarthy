% plot of N brownian motions, compared to   x = sqrt(t)

N = 10;

dt = 0.001;
t = 0:dt:1;

for n=1:N
  x = zeros(size(t));
  for j=2:length(t)
    x(j) = x(j-1) + sqrt(dt)*randn(1);
  end
  plot(t,x)
  hold on
end

plot(t,sqrt(t),'r--','linewidth',6.0)
plot([0 1],[0 0],'k','linewidth',3.0)
plot([0 0],[0 1.5],'k','linewidth',3.0)
hold off
set(gca,'fontsize',36)
xlabel t, ylabel x
axis([0 1 0 1.5]), axis off

print -dpdf brownian.pdf

