# utility module for computing and plotting surface values

import firedrake as fd
from scipy.interpolate import interp1d

# return linear-interpolated z-value along a part of the boundary of the mesh
#   (e.g. surface or base elevation function z = h(x), z = b(x))
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
    xh_halos = fd.Function(P1).interpolate(x).dat.data_with_halos[bc.nodes]  # 1D numpy array
    zh_halos = fd.Function(P1).interpolate(z).dat.data_with_halos[bc.nodes]
    return interp1d(xh_halos,zh_halos,copy=False,
                    bounds_error=False,fill_value='extrapolate')

# return nodal values of spatial coordinate x, and top boundary indicator
# SERIAL
def getxsurface(mesh,top_id):
    P1 = fd.FunctionSpace(mesh, 'CG', 1)
    bc = fd.DirichletBC(P1, 1.0, top_id)
    x,_ = fd.SpatialCoordinate(mesh)
    xh = fd.Function(P1).interpolate(x).dat.data_ro[bc.nodes]
    return (xh, bc)

# return linear-interpolated surface vertical displacement function  r = (delta h)(x)
# SERIAL
def getsurfacevdispfunction(mesh,top_id,r):
    xh,bc = getxsurface(mesh,top_id)
    return interp1d(xh,r.dat.data_ro[bc.nodes],copy=False)

# return linear-interpolated surface velocity functions  u(x), w(x)
# SERIAL
def getsurfacevelocityfunction(mesh,top_id,u):
    xh,bc = getxsurface(mesh,top_id)
    P1V = fd.VectorFunctionSpace(mesh, "CG", 1)
    uP1 = fd.Function(P1V).interpolate(u)
    ufcn = interp1d(xh,uP1.dat.data_ro[bc.nodes,0],copy=False)
    wfcn = interp1d(xh,uP1.dat.data_ro[bc.nodes,1],copy=False)
    return (ufcn,wfcn)

# generate plot of surface values if desired
# SERIAL
def surfaceplot(mesh,u,r,deltat,filename):
    import numpy as np
    import matplotlib.pyplot as plt
    from momentummodel import secpera, dayspera
    from domain import L, bdryids

    x = np.linspace(0.0,L,401)
    hfcn = getboundaryelevation(mesh,bdryids['top'])
    ufcn,wfcn = getsurfacevelocityfunction(mesh,bdryids['top'],u)
    plt.figure(figsize=(6.0,8.0))
    if deltat > 0.0:
        rows = 4
        print('plotting surface values of (h,u,w,h_t) in file %s ...' % filename)
    else:
        rows = 3
        print('plotting surface values of (h,u,w) in file %s ...' % filename)
    plt.subplot(rows,1,1)
    plt.plot(x,hfcn(x),'g',label='surface elevation')
    plt.ylabel('h  [m]')
    plt.legend()

    def removexticks():
        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off

    removexticks()
    plt.subplot(rows,1,2)
    plt.plot(x,secpera*ufcn(x),label='horizontal velocity')
    plt.ylabel('u  [m/a]')
    plt.legend()
    removexticks()
    plt.subplot(rows,1,3)
    plt.plot(x,secpera*wfcn(x),label='vertical velocity')
    plt.ylabel('w  [m/a]')
    plt.legend()
    if deltat > 0.0:
        removexticks()
        rfcn = getsurfacevdispfunction(mesh,bdryids['top'],r)
        plt.subplot(rows,1,4)
        plt.plot(x,rfcn(x)/(deltat/dayspera),'r',label='surface elevation rate (last time step)')
        plt.ylabel('h_t  [m/a]')
        plt.legend()
    plt.xlabel('x  [m]')
    plt.savefig(filename,bbox_inches='tight')

