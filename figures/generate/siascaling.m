% plots versions of Halfar solution to show scaling

J = 500;
[x,y] =  meshgrid(linspace(0,900e3,J+1), [0]);

JJ = 1:5;
secpera = 31556926;
%note t0 = 422.45 a; see halfar.m and Bueler et al (2005)
tt = 10.^(JJ-1) * secpera;  % [1 10 100 1000 10000]  a

figure(1), clf
for j=JJ
  subplot(1,length(JJ),j)
  z = halfar(tt(j),x,y);
  % plot graph and tight box around it
  plot(x/1000,z,'LineWidth',2.0)
  axis([0 900 0 7100])
  set(gca, 'xtick', [0 300 600 900])
  %axis off
end
