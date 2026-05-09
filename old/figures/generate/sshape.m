% plot s-shaped surface not allowed in a shallow ice theory

y = -1.5:.002:1.8;
x = y .* (y.^2 - 1);
bx = linspace(x(1),x(end),20);
plot(x,y,'b','LineWidth',3.0)
hold on
plot(bx,-2.5+0.4*(bx-bx(1))+0.08*randn(size(bx)),'g','LineWidth',3.0)
hold off
axis([x(1) x(end) -2.6 3])
axis off
text(2,-2,'bedrock','FontSize',20)
text(x(end)-2.3,0.5,'ice','FontSize',20)
text(-1,2,'air','FontSize',20)

% to create .eps in Octave
%print -depsc2 sshape

