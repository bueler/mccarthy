#!/usr/bin/env python
from pylab import *

H = 1
u0 = 0.3
z = arange(0.0,H,0.001)
tau = H - z
u = u0 + (H**4.0 - (H-z)**4.0)

zc = arange(0,H+0.1,0.1)
tauc = H - zc
uc = u0 + (H**4.0 - (H-zc)**4.0)

subplot(121)
plot(tau, z, linewidth=3.0), hold(True)
dx = 0.06
dy = 0.02
for k in range(len(zc)-1):
  plot([0,tauc[k]],[zc[k],zc[k]], color='k', linewidth=1.5)
  plot([0,0+dx],[zc[k],zc[k]+dy], color='k', linewidth=1.5)
  plot([tauc[k],tauc[k]-dx],[zc[k],zc[k]-dy], color='k', linewidth=1.5)
plot([-0.2*H,1.2*H],[0,0], color='k')
plot([0,0],[-0.2*H,1.2*H], color='k')
hold(False)
xlabel(r'$\tau_{13}(z)$', fontsize=20.0)
ylabel(r'$z$', fontsize=24.0)
axis([-0.1*H,H,-0.1*H,1.1*H])
setp(gca(), frame_on=False)
setp(gca(), xticks=[], yticks=[0, H], yticklabels = ('0','H'))

subplot(122)
plot(u, z, linewidth=3.0), hold(True)
for k in range(len(zc)):
  plot([0,uc[k]],[zc[k],zc[k]], color='k', linewidth=1.5)
  plot([uc[k],uc[k]-dx],[zc[k],zc[k]+dy], color='k', linewidth=1.5)
  plot([uc[k],uc[k]-dx],[zc[k],zc[k]-dy], color='k', linewidth=1.5)
us = u0+H**4.0 # u scale
plot([-0.2*us,1.1*us],[0,0], color='k')
plot([0,0],[-0.2*H,1.2*H], color='k')
text(u0*0.95,-0.07*H,r'$u_0$',fontsize=16.0)
hold(False)
xlabel(r'$u(z)$', fontsize=20.0)
axis([-0.1*us,1.1*us,-0.1*H,1.1*H])
setp(gca(), frame_on=False)
setp(gca(), yticks=[], xticks=[])


show()

