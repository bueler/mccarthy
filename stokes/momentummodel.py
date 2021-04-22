# module: class MomentumModel

# design principles for this class:
#   1. it owns the velocity, pressure, and mixed spaces
#   2. it does not interact with options, stdout
#   3. it does not know about time stepping or surface kinematical
#   4. it does not interact with files

import sys
import numpy as np
import firedrake as fd

# public data
mixFEchoices = ['P2P1','P3P2','P2P0','CRP0','P1P0']
packagechoices = ['SchurDirect','Direct','SchurGMGSelfp','SchurGMGMass','SchurGMGMassFull']
secpera = 31556926.0    # seconds per year
dayspera = 365.2422     # days per year

# solver packages
# notes on FGMRES:  1. needed if fieldsplit_0,1_ksp_type is *not* preonly
#                   2. always applies PC on right (-s_ksp_pc_side right)
_SchurDirect = {"ksp_type": "fgmres",  # or "gmres" or "minres"
          "pc_type": "fieldsplit",
          "pc_fieldsplit_type": "schur",
          "pc_fieldsplit_schur_factorization_type": "full",  # or "diag"
          "pc_fieldsplit_schur_precondition": "a11",  # the default
          "fieldsplit_0_ksp_type": "preonly",
          "fieldsplit_0_pc_type": "lu",  # uses mumps in parallel
          "fieldsplit_1_ksp_rtol": 1.0e-3,
          "fieldsplit_1_ksp_type": "gmres",
          "fieldsplit_1_pc_type": "none"}
_Direct = {"mat_type": "aij",
          "ksp_type": "preonly",
          "pc_type": "lu",
          "pc_factor_shift_type": "inblocks"}
_SchurGMGSelfp = {"ksp_type": "fgmres",
          "pc_type": "fieldsplit",
          "pc_fieldsplit_type": "schur",
          "pc_fieldsplit_schur_factorization_type": "full",
          "pc_fieldsplit_schur_precondition": "selfp",
          "fieldsplit_0_ksp_type": "preonly",
          "fieldsplit_0_pc_type": "mg",
          "fieldsplit_0_mg_levels_ksp_type": "richardson",
          "fieldsplit_0_mg_levels_pc_type": "sor",  # the default
          "fieldsplit_1_ksp_rtol": 1.0e-3,
          "fieldsplit_1_ksp_type": "gmres",
          "fieldsplit_1_pc_type": "jacobi",
          "fieldsplit_1_pc_jacobi_type": "diagonal"}
# notes on following package:
#   1. right+upper faster than left+lower; not sure why
#   2. -fieldsplit_1_aux_pc_type ilu in serial is implied, but
#      -fieldsplit_1_aux_pc_type bjacobi -fieldsplit_1_aux_sub_pc_type ilu
#      is implied in parallel
#   3. seems to be more robust, but slower, with -fieldsplit_1_aux_pc_type jacobi
#      -fieldsplit_1_aux_pc_jacobi_type rowsum, instead of default ILU
_SchurGMGMass = {"ksp_type": "gmres",
          "ksp_pc_side": "right",
          "pc_type": "fieldsplit",
          "pc_fieldsplit_type": "schur",
          "pc_fieldsplit_schur_factorization_type": "upper",
          "pc_fieldsplit_schur_precondition": "a11",
          "fieldsplit_0_ksp_type": "preonly",
          "fieldsplit_0_pc_type": "mg",
          "fieldsplit_0_mg_levels_ksp_type": "richardson",
          "fieldsplit_0_mg_levels_pc_type": "sor",  # the default
          "fieldsplit_1_ksp_type": "preonly",
          "fieldsplit_1_pc_type": "python",
          "fieldsplit_1_pc_python_type": "momentummodel.Mass"}
          #"fieldsplit_1_aux_pc_type": "jacobi",
          #"fieldsplit_1_aux_pc_jacobi_type": "rowsum"}
_SchurGMGMassFull = {"ksp_type": "fgmres",
          "ksp_pc_side": "right",
          "pc_type": "fieldsplit",
          "pc_fieldsplit_type": "schur",
          "pc_fieldsplit_schur_factorization_type": "full",
          "pc_fieldsplit_schur_precondition": "a11",
          "fieldsplit_0_ksp_type": "preonly",
          "fieldsplit_0_pc_type": "mg",
          "fieldsplit_0_mg_levels_ksp_type": "richardson",
          "fieldsplit_0_mg_levels_pc_type": "sor",  # the default
          "fieldsplit_1_ksp_rtol": 1.0e-3,
          "fieldsplit_1_ksp_type": "gmres",
          "fieldsplit_1_pc_type": "python",
          "fieldsplit_1_pc_python_type": "momentummodel.Mass"}

# match packagechoices order above; default goes first:
packagelist = [_SchurDirect,_Direct,_SchurGMGSelfp,_SchurGMGMass,_SchurGMGMassFull]

# see https://www.firedrakeproject.org/demos/geometric_multigrid.py.html for
# a constant viscosity version of this AuxiliaryOperatorPC approach
class Mass(fd.AuxiliaryOperatorPC):

    def form(self, pc, test, trial):
        # pass in parameters from flow.py
        # FIXME there should be a better way than through the options database?
        opts = fd.PETSc.Options()
        eps = opts.getReal("pcMass_eps")
        Dtyp = opts.getReal("pcMass_Dtyp")
        n_glen = opts.getReal("pcMass_n_glen")
        Bn = opts.getReal("pcMass_Bn")
        more = opts.getReal("pcMass_mass_more_reg")

        # this approach seems to be slightly faster than either keeping u from
        # P2 or using project(u,P1)
        ctx = self.get_appctx(pc)
        u = fd.split(ctx["state"])[0]
        mesh = test.ufl_domain()
        P1 = fd.VectorFunctionSpace(mesh, "CG", 1, dim=2)
        ulow = fd.interpolate(u,P1)
        Du2 = 0.5 * fd.inner(D(ulow), D(ulow)) + (more * eps * Dtyp)**2.0

        # define weak form of viscosity-weighted mass matrix
        mrr = 1.0 - 1.0/n_glen  # positive power if n_glen > 1
        coeff = Du2**(mrr/2.0) / Bn
        a = coeff * fd.inner(test, trial) * fd.dx
        bcs = None
        return (a, bcs)

def D(U):    # strain-rate tensor from velocity U
    return 0.5 * (fd.grad(U) + fd.grad(U).T)

class MomentumModel:

    # physical constants (private)
    _g = 9.81                # m s-2
    _rho = 910.0             # kg m-3
    _A3 = 3.1689e-24         # Pa-3 s-1; EISMINT I value of ice softness
    _B3 = _A3**(-1.0/3.0)    # Pa s(1/3);  ice hardness

    def __init__(self):
        pass

    def set_n_glen(self, n_glen):
        self.n_glen = n_glen

    def get_n_glen(self):
        return self.n_glen

    def set_eps(self, eps):
        self.eps = eps

    def get_eps(self):
        return self.eps

    def set_alpha(self, alpha):
        self.alpha = alpha

    def set_Dtyp_pera(self, Dtyp):
        self.Dtyp = Dtyp / secpera

    def get_Dtyp(self):
        return self.Dtyp

    def set_Hin(self, Hin):
        self.Hin = Hin

    def set_Hout(self, Hout):
        self.Hout = Hout

    # determine B_n so that slab-on-slope solutions have surface velocity
    #   which is n-independent
    def get_Bn(self):
        return (4.0/(self.n_glen+1.0))**(1.0/self.n_glen) \
               * (self._rho*self._g*np.sin(self.alpha)*self.Hin)\
                 **((self.n_glen-3.0)/self.n_glen) \
               * self._B3**(3.0/self.n_glen)

    # compute slab-on-slope inflow velocity
    def _get_uin(self,mesh):
        _,z = fd.SpatialCoordinate(mesh)
        Bn = self.get_Bn()
        C = (2.0 / (self.n_glen + 1.0)) \
            * (self._rho * self._g * np.sin(self.alpha) / Bn)**self.n_glen
        uin = fd.as_vector([C * (self.Hin**(self.n_glen+1.0) \
                                - (self.Hin - z)**(self.n_glen+1.0)),
                            0.0])
        return uin

    def create_mixed_space(self,mesh,mixedtype):
        if   mixedtype == 'P2P1': # Taylor-Hood
            self._V = fd.VectorFunctionSpace(mesh, 'CG', 2)
            self._W = fd.FunctionSpace(mesh, 'CG', 1)
        elif mixedtype == 'P3P2': # Taylor-Hood
            self._V = fd.VectorFunctionSpace(mesh, 'CG', 3)
            self._W = fd.FunctionSpace(mesh, 'CG', 2)
        elif mixedtype == 'P2P0':
            self._V = fd.VectorFunctionSpace(mesh, 'CG', 2)
            self._W = fd.FunctionSpace(mesh, 'DG', 0)
        elif mixedtype == 'CRP0':
            self._V = fd.VectorFunctionSpace(mesh, 'CR', 1)
            self._W = fd.FunctionSpace(mesh, 'DG', 0)
        elif mixedtype == 'P1P0': # interesting but not recommended
            self._V = fd.VectorFunctionSpace(mesh, 'CG', 1)
            self._W = fd.FunctionSpace(mesh, 'DG', 0)
        else:
            print('ERROR: unknown mixed type')
            sys.exit(1)
        self._Z = self._V * self._W

    def solve(self,mesh,bdryids,mixedtype,
              package = 'SchurDirect',
              upold = None, upcoarse = None):
        # define body force and ice hardness
        rhog = self._rho * self._g
        f_body = fd.Constant((rhog * np.sin(self.alpha), - rhog * np.cos(self.alpha)))
        Bn = self.get_Bn()

        # right side outflow nonhomogeneous Neumann is part of weak form;
        #    apply hydrostatic normal force with total force from Hin-height slab
        _,z = fd.SpatialCoordinate(mesh)
        Cout = (self.Hin/self.Hout)**2
        outflow_sigma = fd.as_vector([- Cout * rhog * np.cos(self.alpha) * (self.Hout - z),
                                     Cout * rhog * np.sin(self.alpha) * (self.Hout - z)])

        # create, use old, or prolong coarse solution as initial iterate
        if upold:
            up = upold
        else:
            self.create_mixed_space(mesh,mixedtype)
            up = fd.Function(self._Z)
            if upcoarse:
                fd.prolong(upcoarse,up)
        u,p = fd.split(up)  # get component ufl expressions to define form
        self.u,self.p = up.split()
        self.u.rename('velocity')
        self.p.rename('pressure')

        # define the nonlinear weak form F(u,p;v,q)
        v,q = fd.TestFunctions(self._Z)
        if self.n_glen == 1.0:  # Newtonian ice case
            F = ( fd.inner(Bn * D(u), D(v)) - p * fd.div(v) - fd.div(u) * q \
                  - fd.inner(f_body, v) ) * fd.dx \
                - fd.inner(outflow_sigma, v) * fd.ds(bdryids['outflow'])
        else:                   # (usual) non-Newtonian ice case
            Du2 = 0.5 * fd.inner(D(u), D(u)) + (self.eps * self.Dtyp)**2.0
            rr = 1.0/self.n_glen - 1.0
            F = ( fd.inner(Bn * Du2**(rr/2.0) * D(u), D(v)) \
                  - p * fd.div(v) - fd.div(u) * q - fd.inner(f_body, v) ) * fd.dx \
                - fd.inner(outflow_sigma, v) * fd.ds(bdryids['outflow'])

        # Dirichlet boundary conditions
        noslip = fd.Constant((0.0, 0.0))
        inflow_u = self._get_uin(mesh)
        bcs = [ fd.DirichletBC(self._Z.sub(0), noslip, bdryids['base']),
                fd.DirichletBC(self._Z.sub(0), inflow_u, bdryids['inflow']) ]

        # solve
        params = packagelist[packagechoices.index(package)]
        # note: Sometime in 2020, Firedrake decided to override the default
        #       PETSc linesearch of 'bt' and replace it with 'basic', i.e. NO
        #       linesearch.  We need linesearch for n = 3 etc.
        params['snes_linesearch_type'] = 'bt'
        fd.solve(F == 0, up, bcs=bcs, options_prefix='s',
                 solver_parameters=params)

        return up

    # statistics about the solution:
    #   umagav  = average velocity magnitude
    #   umagmax = maximum velocity magnitude
    #   pav     = average pressure
    #   pmax    = maximum pressure
    def solutionstats(self,mesh):
        P1 = fd.FunctionSpace(mesh, 'CG', 1)
        one = fd.Constant(1.0, domain=mesh)
        area = fd.assemble(fd.dot(one,one) * fd.dx)
        pav = fd.assemble(fd.sqrt(fd.dot(self.p, self.p)) * fd.dx) / area
        with self.p.dat.vec_ro as vp:
            pmax = vp.max()[1]
        umagav = fd.assemble(fd.sqrt(fd.dot(self.u, self.u)) * fd.dx) / area
        umag = fd.interpolate(fd.sqrt(fd.dot(self.u,self.u)),P1)
        with umag.dat.vec_ro as vumag:
            umagmax = vumag.max()[1]
        return umagav,umagmax,pav,pmax

    # generate regularized effective viscosity from the solution:
    #   nu = (1/2) B_n X^((1/n)-1)
    # where X = sqrt(|Du|^2 + eps^2 Dtyp^2)
    def effectiveviscosity(self,mesh):
        P1 = fd.FunctionSpace(mesh, 'CG', 1)
        Du2 = 0.5 * fd.inner(D(self.u), D(self.u)) + (self.eps * self.Dtyp)**2.0
        rr = 1.0/self.n_glen - 1.0
        Bn = self.get_Bn()
        nu = fd.interpolate(0.5 * Bn * Du2**(rr/2.0), P1)
        nu.rename('effective viscosity')
        return nu

    # numerical errors relative to slab-on-slope solution:
    #   uerrmax = maximum magnitude of velocity error
    #   perrmax = maximum of pressure error
    def numerical_errors_slab(self,mesh):
        P1 = fd.FunctionSpace(mesh, 'CG', 1)
        up_exact = fd.Function(self._Z)
        u_exact,p_exact = up_exact.split()
        inflow_u = self._get_uin(mesh)
        u_exact.interpolate(inflow_u)
        _,z = fd.SpatialCoordinate(mesh)
        p_exact.interpolate(self._rho * self._g * np.cos(self.alpha) \
                            * (self.Hin - z))
        uerr = fd.interpolate(fd.sqrt(fd.dot(u_exact-self.u,u_exact-self.u)),P1)
        perr = fd.interpolate(fd.sqrt(fd.dot(p_exact-self.p,p_exact-self.p)),self._W)
        with uerr.dat.vec_ro as vuerr:
            uerrmax = vuerr.max()[1]
        with perr.dat.vec_ro as vperr:
            perrmax = vperr.max()[1]
        return uerrmax,perrmax

