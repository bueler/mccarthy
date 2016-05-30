#!/usr/bin/env python
from pylab import *

x = arange(0.0,5.0,0.01);
H0 = 1000.0
H = H0 + 0.1 * H0 * sin(x**1.5)
plot(x, H, color='k', linewidth=2.0), hold(True)
plot([-0.5,5.5],[0.0,0.0], color='k')
plot([-0.5,],[-0.1], color='k', linewidth=0.001)
plot([0.0,5.0],[0.0,0.0], color='k', linewidth=2.0)
plot([0.0,0.0],[0.0,H0], color='k', linewidth=2.0)
plot([5.0,5.0],[0.0,H[-1]], color='k', linewidth=2.0)

annotate(r'$\bar U_1$', size=24.0,
         xy=(0.5, H0/3.0), xytext=(-0.6, H0/3.0-0.05*H0), arrowprops=dict(arrowstyle="->") )
annotate(r'$\bar U_2$', size=24.0,
         xy=(5.5, H0/3.0), xytext=(4.0, H0/3.0-0.05*H0), arrowprops=dict(arrowstyle="->") )

annotate(r'$M(x)$', size=24.0,
         xy=(2.5, H[200]), xytext=(2.1, 1.15 * H[200]), arrowprops=dict(arrowstyle="->") )

text(2.5,H0*0.65,r'$A$', size=24.0)
text(-0.5,0.7*H0,r'$H_1$', size=24.0)
text(5.1,0.7*H[-1],r'$H_2$', size=24.0)

text(-0.1,-0.15*H0,r'$x_1$', size=24.0)
text(4.9,-0.15*H0,r'$x_2$', size=24.0)

hold(False)

axis('off')

show()

