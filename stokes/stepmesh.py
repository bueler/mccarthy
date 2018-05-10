# module for glacier mesh metadata

from firedrake import *

# numbering of parts of boundary *must match generation script genstepmesh.py*
bdryids = {'outflow' : 41,
           'top'     : 42,
           'inflow'  : 43,
           'base'    : 44}

# extract mesh geometry making these definitions (tolerance=1cm):
#   L       = total length of domain     = (max x-coordinate)
#   bs      = height of bed step         = (min z-coordinate at x=0)
#   hsurfin = height of surface at input = (max z-coordinate at x=0)
#   Hout    = ice thickness at output    = (max z-coordinate at x=L)
# in parallel no process owns the whole mesh so MPI_Allreduce() is needed
from mpi4py import MPI
def getmeshdims(mesh,tol=0.01):
    dims = {}
    xa = mesh.coordinates.dat.data_ro[:,0]  # .data_ro acts like VecGetArrayRead
    za = mesh.coordinates.dat.data_ro[:,1]
    loc_L = max(xa)
    loc_bs = 9.9999e99
    loc_hsurfin = 0.0
    if any(xa < tol):              # some processes may have no such points
        loc_bs = min(za[xa<tol])
        loc_hsurfin = max(za[xa<tol])
    L = mesh.comm.allreduce(loc_L, op=MPI.MAX)
    dims['L'] = L
    dims['bs'] = mesh.comm.allreduce(loc_bs, op=MPI.MIN)
    dims['hsurfin'] = mesh.comm.allreduce(loc_hsurfin, op=MPI.MAX)
    # note that determining Hout requires already having L
    loc_Hout = 0.0
    if any(xa > L-tol):
        loc_Hout = max(za[xa > L-tol])
    dims['Hout'] = mesh.comm.allreduce(loc_Hout, op=MPI.MAX)
    dims['isslab'] = False
    Hin = dims['hsurfin'] - dims['bs']
    if abs(dims['bs']) < tol and abs(Hin - dims['Hout']) < tol:
        dims['isslab']= True
    return dims

