# module provides writeout()

import matplotlib.pyplot as plt

_SHOW = False
def writeout(outname):
    if _SHOW:
        plt.show()
    else:
        print('writing file ' + outname)
        plt.savefig(outname, bbox_inches='tight', transparent=True)
