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

# return nodal values of spatial coordinate z, and top boundary indicator
# SERIAL
def getzsurface(mesh,top_id):
    P1 = fd.FunctionSpace(mesh, 'CG', 1)
    bc = fd.DirichletBC(P1, 1.0, top_id)
    _,z = fd.SpatialCoordinate(mesh)
    zs = fd.Function(P1).interpolate(z).dat.data_ro[bc.nodes]
    return (zs, bc)

# return linear-interpolated surface vertical displacement function  r = (delta h)(x)
# SERIAL
def getsurfacevdispfunction(mesh,top_id,r):
    xs,bc = getxsurface(mesh,top_id)
    return interp1d(xs,r.dat.data_ro[bc.nodes],copy=False)

# return linear-interpolated surface velocity functions  u(x), w(x)
# SERIAL
def getsurfacevelocityfunction(mesh,top_id,u):
    xs,bc = getxsurface(mesh,top_id)
    P1V = fd.VectorFunctionSpace(mesh, 'CG', 1)
    uP1 = fd.Function(P1V).interpolate(u)
    ufcn = interp1d(xs,uP1.dat.data_ro[bc.nodes,0],copy=False)
    wfcn = interp1d(xs,uP1.dat.data_ro[bc.nodes,1],copy=False)
    return (ufcn,wfcn)

# FIXME this is a kludge version which averages u,w to same staggered locations as dsdx
# return values of  - vecu|_s . n_s = + u_s s_x - w_s  as 1D numpy array
# with one value for each surface mesh edge; note n_s = <-s_x,1>
# SERIAL
def getbalancemotion(mesh,top_id,u):
    import numpy as np
    xs,bc = getxsurface(mesh,top_id)
    if np.min(np.abs(np.diff(xs))) < 0.001:
        print('WARNING: mesh has |delta x| < 1 mm at surface')
    zs,_ = getzsurface(mesh,top_id)
    P1V = fd.VectorFunctionSpace(mesh, 'CG', 1)
    uP1 = fd.Function(P1V).interpolate(u)
    us = uP1.dat.data_ro[bc.nodes,0]
    ws = uP1.dat.data_ro[bc.nodes,1]
    xsav = (xs[:-1] + xs[1:]) / 2.0
    dsdx = np.diff(zs) / np.diff(xs)
    usav = (us[:-1] + us[1:]) / 2.0
    wsav = (ws[:-1] + ws[1:]) / 2.0
    return xsav, usav * dsdx - wsav

def removexticks():
    import matplotlib.pyplot as plt
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

# generate plot of surface values if desired
# SERIAL
def surfaceplot(mesh,u,r,deltat,filename):
    import numpy as np
    import matplotlib.pyplot as plt
    from momentummodel import secpera, dayspera
    from domain import getdomaindims, bdryids

    L,_,_,_,_ = getdomaindims(mesh)
    x = np.linspace(0.0,L,401)
    sfcn = getboundaryelevation(mesh,bdryids['top'])
    ufcn,wfcn = getsurfacevelocityfunction(mesh,bdryids['top'],u)
    plt.figure(figsize=(6.0,8.0))
    if deltat > 0.0:
        rows = 5
        print('plotting surface values of (s,u,w,Phi,s_t) in file %s ...' % filename)
    else:
        rows = 4
        print('plotting surface values of (s,u,w,Phi) in file %s ...' % filename)
    plt.subplot(rows,1,1)
    plt.plot(x,sfcn(x),'g',label='surface elevation')
    plt.ylabel('s  [m]')
    plt.legend()
    removexticks()
    plt.subplot(rows,1,2)
    plt.plot(x,secpera*ufcn(x),label='horizontal velocity')
    plt.ylabel('u  [m/a]')
    plt.legend()
    removexticks()
    plt.subplot(rows,1,3)
    plt.plot(x,secpera*wfcn(x),label='vertical velocity')
    plt.plot(x,np.zeros(np.shape(x)),'k-',lw=0.5)
    plt.ylabel('w  [m/a]')
    plt.legend()
    xs, Phi = getbalancemotion(mesh,bdryids['top'],u)
    removexticks()
    plt.subplot(rows,1,4)
    plt.plot(xs,secpera*Phi,'.',label='balance motion')
    plt.plot(xs,np.zeros(np.shape(xs)),'k-',lw=0.5)
    plt.ylabel('Phi  [m/a]')
    plt.legend()
    if deltat > 0.0:
        removexticks()
        rfcn = getsurfacevdispfunction(mesh,bdryids['top'],r)
        plt.subplot(rows,1,5)
        plt.plot(x,rfcn(x)/(deltat/dayspera),'r',label='surface elevation rate (last time step)')
        plt.ylabel('s_t  [m/a]')
        plt.legend()
    plt.xlabel('x  [m]')
    plt.savefig(filename,bbox_inches='tight')

def surfaceplotgreen(mesh,ugreen,greenx,filename):
    import numpy as np
    import matplotlib.pyplot as plt
    from momentummodel import secpera, dayspera
    from domain import getdomaindims, bdryids

    L,_,_,_,_ = getdomaindims(mesh)
    dxdrop = 20.0
    assert (greenx > dxdrop)
    assert (greenx < L - dxdrop)
    xleft = np.linspace(0.0,greenx - dxdrop,401)
    xright = np.linspace(greenx + dxdrop,L,401)
    sfcn = getboundaryelevation(mesh,bdryids['top'])
    ufcn,wfcn = getsurfacevelocityfunction(mesh,bdryids['top'],ugreen)
    plt.figure(figsize=(6.0,5.0))
    print("plotting Green's function surface values of (u,w) in file %s ..." % filename)
    plt.subplot(2,1,1)
    plt.plot(greenx,0.0,'k.',ms=8.0)
    plt.plot([0.0,L],[0.0,0.0],'k--',lw=0.5)
    plt.plot(xleft,secpera*ufcn(xleft),'g',label='horizontal velocity')
    plt.plot(xright,secpera*ufcn(xright),'g')
    plt.ylabel('u  [m/a]')
    plt.legend()
    removexticks()
    plt.subplot(2,1,2)
    plt.plot(greenx,0.0,'k.',ms=8.0)
    plt.plot([0.0,L],[0.0,0.0],'k--',lw=0.5)
    plt.plot(xleft,secpera*wfcn(xleft),'g',label='vertical velocity')
    plt.plot(xright,secpera*wfcn(xright),'g')
    plt.ylabel('w  [m/a]')
    plt.legend()
    plt.xlabel('x  [m]')
    plt.savefig(filename,bbox_inches='tight')
