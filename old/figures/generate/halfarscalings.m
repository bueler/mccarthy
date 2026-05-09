% plot of powers of t arising in Halfar solution

L = 3;
dt = 0.002;
t = dt:dt:L;

plot(t,t.^(-1/9),'r','linewidth',3.0)
hold on
plot(t,t.^(-1/18),'b','linewidth',3.0)
plot([0 L],[1 1],'k--','linewidth',1.5)
grid on
set(gca,'fontsize',20)
xlabel t

axis([0 L 0 2])

print -dpdf halfarscalings.pdf

