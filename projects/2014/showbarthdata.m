load Bartholomaus2011_Fig7.txt

t = Bartholomaus2011_Fig7(:,1);
Q_in = Bartholomaus2011_Fig7(:,2);
Q_out = Bartholomaus2011_Fig7(:,3);
U_sliding = Bartholomaus2011_Fig7(:,4);

fprintf('data has %d records, from day %.2f to day %.2f ...\n',length(t),min(t),max(t))

subplot(2,1,1)
plot(t,Q_in,'b',t,Q_out,'r')
legend('Q Input','Q Output')
ylabel('Water flux  (m^3 s^{-1})')
axis([130 230 0 1000])

subplot(2,1,2)
plot(t,U_sliding,'k')
xlabel('Day of year 2006')
ylabel('u_b  (m d^{-1})')
axis([130 230 0 2])

