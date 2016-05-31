#!/usr/bin/env python
from pylab import *

L = 1.0
T = 1.0
t = linspace(0.0,T,11)
x = linspace(0.0,L,7)

figure(figsize=(4.0,5.0))
hold(True)

# draw axes
plot([-0.04,1.1*L],[0.0,0.0], color='k')
plot([0.0,0.0],[-0.03,1.1*T], color='k')
text(1.15,0.0,r'$x$', fontsize=16.0)
text(0.0,1.15,r'$t$', fontsize=18.0)
plot([L,L],[-0.03,0.0], color='k')
plot([-0.04,0.0], [L,L], color='k')
#text(-0.04,-0.15,r'$0$', fontsize=20.0)
#text(0.95*L,-0.15,r'$L$', fontsize=20.0)
#text(-0.15,-0.04,r'$0$', fontsize=20.0)
#text(-0.15,0.97*T,r'$T$', fontsize=20.0)

# draw grid
for xj in x:
    plot([xj,xj],[0.0,T], color='k', linewidth=1.5)  # vertical lines
for tn in t:
    plot([0.0,L],[tn,tn], color='k', linewidth=1.5)  # horizontal lines

# label generic point
plot(x[4],t[7],'ko')
text(0.93*x[4],-0.1,r'$x_j$', fontsize=20.0)
text(-0.15,0.95*t[7],r'$t_n$', fontsize=20.0)

axis('off')

#show()
savefig('timespacegrid.pdf')

