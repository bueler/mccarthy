import numpy as np
import matplotlib.pyplot as plt
from writeout import *

x = np.arange(0.0,5.0,0.01);
H0 = 1000.0
H = H0 + 0.1 * H0 * np.sin(x**1.5)
plt.plot(x, H, color='k', linewidth=2.0)
plt.plot([-0.5,5.5],[0.0,0.0], color='k')
plt.plot([-0.5,],[-0.1], color='k', linewidth=0.001)
plt.plot([0.0,5.0],[0.0,0.0], color='k', linewidth=2.0)
plt.plot([0.0,0.0],[0.0,H0], color='k', linewidth=2.0)
plt.plot([5.0,5.0],[0.0,H[-1]], color='k', linewidth=2.0)

plt.annotate("", xy=(0.8, H0/3.0), xytext=(-0.9, H0/3.0), arrowprops=dict(arrowstyle="->") )
plt.annotate("", xy=(5.5, H0/3.0), xytext=(4.5, H0/3.0), arrowprops=dict(arrowstyle="->") )

#plt.annotate("", xy=(2.5, H[200]), xytext=(2.1, 1.15 * H[200]), arrowprops=dict(arrowstyle="->"))

#plt.text(2.5,H0*0.65,'(area)', size=24.0)
plt.text(2.3,H0*1.05,r"$a$", size=24.0)

plt.text(-0.6,0.6*H0,r'$H_1$', size=24.0)
plt.text(5.1,0.6*H0,r'$H_2$', size=24.0)

plt.text(-0.6,0.15*H0,r'$\bar U_1$', size=24.0)
plt.text(5.1,0.15*H0,r'$\bar U_2$', size=24.0)

plt.text(-0.1,-0.15*H0,r'$x_1$', size=24.0)
plt.text(4.9,-0.15*H0,r'$x_2$', size=24.0)

plt.axis('off')
#plt.show()
writeout("slabmasscontfig.pdf")
