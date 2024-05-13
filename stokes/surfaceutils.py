# utility module for computing and plotting surface values

import firedrake as fd
from scipy.interpolate import interp1d

# return linear-interpolated z-value along a part of the boundary of the mesh
#   (e.g. surface or base elevation function z = s(x), z = b(x))
# PARALLEL
# notes: 1)  bc.nodes  gives indices to mesh; is a 1D dtype=int32 numpy array
#        2)  bc.nodes  includes the halo nodes so need .dat.data_with_halos
#        3)  f = Function(P1);  bc.apply(f)  would give indicator function
#        4)  default for interp1d() is 'linear', which is what we want
#        5)  by turning off bounds_error, and allowing extrapolation,
#            there is parallel functionality for this function
def getboundaryelevation(mesh,bdrypartid):
    P1 = fd.FunctionSpace(mesh, 'CG', 1)
    bc = fd.DirichletBC(P1, 1.0, bdrypartid)
    x,z = fd.SpatialCoordinate(mesh)
    xs_halos = fd.Function(P1).interpolate(x).dat.data_with_halos[bc.nodes]  # 1D numpy array
    zs_halos = fd.Function(P1).interpolate(z).dat.data_with_halos[bc.nodes]
    return interp1d(xs_halos,zs_halos,copy=False,
                    bounds_error=False,fill_value='extrapolate')

# return nodal values of spatial coordinate x, and top boundary indicator
# SERIAL
def getxsurface(mesh,top_id):
    P1 = fd.FunctionSpace(mesh, 'CG', 1)
    bc = fd.DirichletBC(P1, 1.0, top_id)
    x,_ = fd.SpatialCoordinate(mesh)
    xs = fd.Function(P1).interpolate(x).dat.data_ro[bc.nodes]
    return (xs, bc)

# return linear-interpolated surface velocity functions  u(x), w(x)
# SERIAL
def getsurfacevelocityfunction(mesh,top_id,u):
    xs,bc = getxsurface(mesh,top_id)
    P1V = fd.VectorFunctionSpace(mesh, 'CG', 1)
    uP1 = fd.Function(P1V).interpolate(u)
    ufcn = interp1d(xs,uP1.dat.data_ro[bc.nodes,0],copy=False)
    wfcn = interp1d(xs,uP1.dat.data_ro[bc.nodes,1],copy=False)
    return (ufcn, wfcn)

def removexticks():
    import matplotlib.pyplot as plt
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

# generate plot of surface velocity values; SERIAL
def surfaceplot(mesh, u, filename):
    import numpy as np
    import matplotlib.pyplot as plt
    from momentummodel import secpera, dayspera
    from domain import bdryids

    xs, _ = getxsurface(mesh, bdryids['top'])
    L = max(xs)
    x = np.linspace(0.0,L,401)
    sfcn = getboundaryelevation(mesh, bdryids['top'])
    ufcn, wfcn = getsurfacevelocityfunction(mesh, bdryids['top'], u)
    plt.figure(figsize=(6.0,8.0))
    rows = 2
    print('plotting surface values of (u,w) in image file %s ...' % filename)
    plt.subplot(rows,1,1)
    plt.plot(x,secpera*ufcn(x),label='horizontal velocity')
    plt.ylabel('u  [m/a]')
    plt.legend()
    removexticks()
    plt.subplot(rows,1,2)
    plt.plot(x,secpera*wfcn(x),label='vertical velocity')
    plt.plot(x,np.zeros(np.shape(x)),'k-',lw=0.5)
    plt.ylabel('w  [m/a]')
    plt.legend()
    plt.xlabel('x  [m]')
    plt.savefig(filename,bbox_inches='tight')
