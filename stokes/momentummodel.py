# module defining class MomentumModel

# design principles for this class:
#  -1. FIXME THIS IS BAD IDEA?  it owns the mesh
#   0. it owns the velocity,pressure spaces
#   1. it does not interact with options, stdout
#   2. it does not know about time stepping or surface kinematical
#   3. FIXME for now it calls existing functionality in physics.py, meshactions.py
#   4. it does not interact with files

from firedrake import *
import sys
from physics import stokessolve, solutionstats, numericalerrors_slab

mixFEchoices = ['P2P1','P3P2','P2P0','CRP0','P1P0']

class MomentumModel:

    _secpera = 31556926.0       # seconds per year
    #FIXME following are in physics.py for now
    #g = 9.81                   # m s-2
    #rho = 910.0                # kg m-3
    #A3 = 3.1689e-24            # Pa-3 s-1; EISMINT I value of ice softness
    #B3 = A3**(-1.0/3.0)        # Pa s(1/3);  ice hardness

    # store mesh and define mixed finite element space
    def __init__(self, mesh, bdryids, mixedtype):
        self.mesh = mesh
        self.bdryids = bdryids
        if   mixedtype == 'P2P1': # Taylor-Hood
            self.V = VectorFunctionSpace(self.mesh, 'CG', 2)
            self.W = FunctionSpace(self.mesh, 'CG', 1)
        elif mixedtype == 'P3P2': # Taylor-Hood
            self.V = VectorFunctionSpace(self.mesh, 'CG', 3)
            self.W = FunctionSpace(self.mesh, 'CG', 2)
        elif mixedtype == 'P2P0':
            self.V = VectorFunctionSpace(self.mesh, 'CG', 2)
            self.W = FunctionSpace(self.mesh, 'DG', 0)
        elif mixedtype == 'CRP0':
            self.V = VectorFunctionSpace(self.mesh, 'CR', 1)
            self.W = FunctionSpace(self.mesh, 'DG', 0)
        elif mixedtype == 'P1P0': # interesting but not recommended
            self.V = VectorFunctionSpace(self.mesh, 'CG', 1)
            self.W = FunctionSpace(self.mesh, 'DG', 0)
        else:
            print('ERROR: unknown mixed type')
            sys.exit(1)
        self.Z = self.V * self.W

    def secpera(self):
        return self._secpera

    # FIXME self.mesh should be private!

    def set_mesh_coordinates(self,f):
        self.mesh.coordinates.assign(f)
        bmin_initial = 0.0  # FIXME
        # mesh z values below bed is extreme instability
        return any(f.dat.data_ro[:,1] < bmin_initial - 1.0)

    # FIXME self.Z should be private!

    def set_nglen(self, n_glen):
        self.n_glen = n_glen

    def set_eps(self, eps):
        self.eps = eps

    def set_alpha(self, alpha):
        self.alpha = alpha

    def set_Dtyp_pera(self, Dtyp):
        self.Dtyp = Dtyp / self._secpera

    def set_Hin(self, Hin):
        self.Hin = Hin

    def set_Hout(self, Hout):
        self.Hout = Hout

    def solve(self):  # FIXME just calls physics
        self.up = stokessolve(Function(self.Z),
                              self.mesh,self.bdryids,self.Z,
                              Hin = self.Hin,
                              Hout = self.Hout,
                              n_glen = self.n_glen,
                              alpha = self.alpha,
                              eps = self.eps,
                              Dtyp = self.Dtyp)

    def getsolution(self):
        u,p = self.up.split()
        u.rename('velocity')
        p.rename('pressure')
        return u,p

    def solutionstats(self):  # FIXME just calls physics
        u,p = self.getsolution()
        umagav,umagmax,pav,pmax = solutionstats(u,p,self.mesh)
        return umagav,umagmax,pav,pmax

    def numerical_errors_slab(self):  # FIXME just calls physics
        u,p = self.getsolution()
        uerrmax,perrmax = numericalerrors_slab(u,p,self.mesh,self.V,self.W,
                                               self.Hin,self.n_glen,self.alpha)
        return uerrmax,perrmax


