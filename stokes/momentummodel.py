# module: class MomentumModel

# design principles for this class:
#   1. it owns the velocity, pressure, and mixed spaces
#   2. it does not interact with options or stdout
#   3. it does not know about time stepping or surface kinematical equation
#   4. it does not interact with files

import sys
import numpy as np
import firedrake as fd
from firedrake.petsc import PETSc, OptionsManager

# public data
secpera = 31556926.0    # seconds per year
dayspera = 365.2422     # days per year

# solver packages
# notes on FGMRES:  1. needed if fieldsplit_0,1_ksp_type is *not* preonly
#                   2. always applies PC on right (-s_ksp_pc_side right)
_SchurDirect = {"ksp_type": "fgmres",
          "pc_type": "fieldsplit",
          "pc_fieldsplit_type": "schur",
          "pc_fieldsplit_schur_factorization_type": "full",  # or "diag"
          "pc_fieldsplit_schur_precondition": "a11",  # the default
          "fieldsplit_0_ksp_type": "preonly",
          "fieldsplit_0_pc_type": "lu",
          "fieldsplit_0_pc_factor_mat_solver_type": "mumps",
          "fieldsplit_1_ksp_rtol": 1.0e-3,
          "fieldsplit_1_ksp_type": "gmres",
          "fieldsplit_1_pc_type": "none"}
_Direct = {"mat_type": "aij",
          "ksp_type": "preonly",
          "pc_type": "lu",
          "pc_factor_shift_type": "inblocks",
          "pc_factor_mat_solver_type": "mumps"}

packagedict = {"SchurDirect": _SchurDirect,
               "Direct":      _Direct}

def D(U):    # strain-rate tensor from velocity U
    return 0.5 * (fd.grad(U) + fd.grad(U).T)

class MomentumModel(OptionsManager):

    def __init__(self, **kwargs):
        # physical constants (private)
        self._g = 9.81                # m s-2
        self._rho = 910.0             # kg m-3
        self._A3 = 3.1689e-24         # Pa-3 s-1; EISMINT I value of ice softness
        self._B3 = self._A3**(-1.0/3.0)    # Pa s(1/3);  ice hardness
        self.n_glen = 3.0
        self.eps = kwargs.pop("eps", 0.01)
        self.alpha = kwargs.pop("alpha", 0.1)
        self.Dtyp = kwargs.pop("Dtyp_pera", 2.0) / secpera
        self.Hin = kwargs.pop("Hin", 400.0)
        self.Hout = kwargs.pop("Hout", 400.0)

    # compute slab-on-slope inflow velocity; note alpha = 0 ==> uin = 0
    def _get_uin(self,mesh):
        _,z = fd.SpatialCoordinate(mesh)
        C = (2.0 / (self.n_glen + 1.0)) \
            * (self._rho * self._g * np.sin(self.alpha) / self._B3)**self.n_glen
        uin = fd.as_vector([C * (self.Hin**(self.n_glen+1.0) \
                                - (self.Hin - z)**(self.n_glen+1.0)),
                            0.0])
        return uin

    def create_mixed_space(self, mesh):
        self._V = fd.VectorFunctionSpace(mesh, 'CG', 2)
        self._W = fd.FunctionSpace(mesh, 'CG', 1)
        self._Z = self._V * self._W

    def solve(self, mesh, bdryids, package = 'SchurDirect',
              upold = None, upcoarse = None):
        # define body force and ice hardness
        rhog = self._rho * self._g
        f_body = fd.Constant((rhog * np.sin(self.alpha), - rhog * np.cos(self.alpha)))

        if self.Hout >= 1.0:
            # if there is an outflow boundary then it is Neumann
            # and thus part of the weak form; here we apply nonhomogeneous:
            # hydrostatic normal force equivalent to height=Hin slab,
            # as though there was a down-stream glacier continuation
            _,z = fd.SpatialCoordinate(mesh)
            Cout = (self.Hin / self.Hout)**2
            outflow_sigma = fd.as_vector( \
                [- Cout * rhog * np.cos(self.alpha) * (self.Hout - z),
                Cout * rhog * np.sin(self.alpha) * (self.Hout - z)])

        # create, use old, or prolong coarse solution as initial iterate
        if upold:
            up = upold
        else:
            self.create_mixed_space(mesh)
            up = fd.Function(self._Z)
            if upcoarse:
                fd.prolong(upcoarse,up)
        u, p = fd.split(up)  # get component ufl expressions to define form
        self.u, self.p = up.subfunctions
        self.u.rename('velocity')
        self.p.rename('pressure')

        # define the nonlinear weak form F(u,p;v,q)
        v,q = fd.TestFunctions(self._Z)
        if self.n_glen == 1.0:  # Newtonian ice case
            F = ( fd.inner(self._B3 * D(u), D(v)) - p * fd.div(v) - fd.div(u) * q \
                  - fd.inner(f_body, v) ) * fd.dx
        else:                   # (usual) non-Newtonian ice case
            Du2 = 0.5 * fd.inner(D(u), D(u)) + (self.eps * self.Dtyp)**2.0
            rr = 1.0/self.n_glen - 1.0
            F = ( fd.inner(self._B3 * Du2**(rr/2.0) * D(u), D(v)) \
                  - p * fd.div(v) - fd.div(u) * q - fd.inner(f_body, v) ) * fd.dx
        if self.Hout >= 1.0:
            F = F - fd.inner(outflow_sigma, v) * fd.ds(bdryids['outflow'])

        # Dirichlet boundary conditions
        noslip = fd.Constant((0.0, 0.0))
        inflow_u = self._get_uin(mesh)
        bcs = [ fd.DirichletBC(self._Z.sub(0), noslip, bdryids['base']),
                fd.DirichletBC(self._Z.sub(0), inflow_u, bdryids['inflow']) ]

        # solve
        params = packagedict[package]
        # note: the Firedrake default is 'basic', i.e. no linesearch
        #       we need 'bt' linesearch for n = 3
        params['snes_linesearch_type'] = 'bt'
        fd.solve(F == 0, up, bcs=bcs, options_prefix='s',
                 solver_parameters=params)

        return up

    # statistics about the solution:
    #   umagav  = average velocity magnitude
    #   umagmax = maximum velocity magnitude
    #   pav     = average pressure
    #   pmax    = maximum pressure
    def solutionstats(self, mesh):
        P1 = fd.FunctionSpace(mesh, 'CG', 1)
        R = fd.FunctionSpace(mesh, 'R', 0)
        one = fd.Function(R).assign(1.0)
        area = fd.assemble(one * fd.dx)
        pav = fd.assemble(fd.sqrt(fd.dot(self.p, self.p)) * fd.dx) / area
        with self.p.dat.vec_ro as vp:
            pmax = vp.max()[1]
        umagav = fd.assemble(fd.sqrt(fd.dot(self.u, self.u)) * fd.dx) / area
        umag = fd.Function(P1).interpolate(fd.sqrt(fd.dot(self.u,self.u)))
        with umag.dat.vec_ro as vumag:
            umagmax = vumag.max()[1]
        return umagav,umagmax,pav,pmax

    # generate regularized effective viscosity from the solution:
    #   nu = (1/2) B_n X^((1/n)-1)
    # where X = sqrt(|Du|^2 + eps^2 Dtyp^2)
    def effectiveviscosity(self, mesh):
        P1 = fd.FunctionSpace(mesh, 'CG', 1)
        Du2 = 0.5 * fd.inner(D(self.u), D(self.u)) + (self.eps * self.Dtyp)**2.0
        rr = 1.0 / self.n_glen - 1.0
        nu = fd.Function(P1).interpolate(0.5 * self._B3 * Du2**(rr/2.0))
        nu.rename('effective viscosity')
        return nu

    # numerical errors relative to slab-on-slope solution:
    #   uerrmax = maximum magnitude of velocity error
    #   perrmax = maximum of pressure error
    def numerical_errors_slab(self,mesh):
        P1 = fd.FunctionSpace(mesh, 'CG', 1)
        up_exact = fd.Function(self._Z)
        u_exact, p_exact = up_exact.subfunctions
        inflow_u = self._get_uin(mesh)
        u_exact.interpolate(inflow_u)
        _,z = fd.SpatialCoordinate(mesh)
        p_exact.interpolate(self._rho * self._g * np.cos(self.alpha) \
                            * (self.Hin - z))
        uerr = fd.Function(P1).interpolate(fd.sqrt(fd.dot(u_exact-self.u,u_exact-self.u)))
        perr = fd.Function(self._W).interpolate(fd.sqrt(fd.dot(p_exact-self.p,p_exact-self.p)))
        with uerr.dat.vec_ro as vuerr:
            uerrmax = vuerr.max()[1]
        with perr.dat.vec_ro as vperr:
            perrmax = vperr.max()[1]
        return uerrmax,perrmax
