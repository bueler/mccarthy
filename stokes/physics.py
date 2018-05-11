# module for glacier stokes solver physics

from firedrake import *

secpera = 31556926.0       # seconds per year

# glacier physical constants
g = 9.81                   # m s-2
rho = 910.0                # kg m-3
A3 = 3.1689e-24            # Pa-3 s-1; EISMINT I value of ice softness
B3 = A3**(-1.0/3.0)        # Pa s(1/3);  ice hardness

# determine B_n so that slab-on-slope solutions give surface velocity that is n-independent
def getBn(n_glen,alpha,Hin):
    Bn = (4.0/(n_glen+1.0))**(1.0/n_glen) \
         * (rho*g*sin(alpha)*Hin)**((n_glen-3.0)/n_glen) * B3**(3.0/n_glen)
    return Bn

# define the strain-rate tensor from the velocity U
def D(U):
    return 0.5 * (grad(U) + grad(U).T)

# compute slab-on-slope inflow velocity
def getinflow(mesh,hsurfin,Hin,n_glen,alpha):
    _,z = SpatialCoordinate(mesh)
    Bn = getBn(n_glen,alpha,Hin)
    C = (2.0 / (n_glen + 1.0)) * (rho * g * sin(alpha) / Bn)**n_glen
    u = as_vector([C * (Hin**(n_glen+1.0) - (hsurfin - z)**(n_glen+1.0)), 0.0])
    return u

def stokessolve(mesh,bdryids,Z,hsurfin,Hin,Hout,n_glen,alpha,eps,Dtyp):
    # define body force
    f_body = Constant((g * rho * sin(alpha), - g * rho * cos(alpha)))

    # define ice hardness
    Bn = getBn(n_glen,alpha,Hin)

    # right side outflow nonhomogeneous Neumann is part of weak form:
    #    apply hydrostatic normal force; 
    x,z = SpatialCoordinate(mesh)
    outflow_sigma = as_vector([- rho * g * cos(alpha) * (Hout - z),
                               rho * g * sin(alpha) * (Hout - z)])

    # define the nonlinear weak form F(u,p;v,q)
    up = Function(Z)       # *not* TrialFunctions(Z)
    u,p = split(up)        # up.split() not equivalent here?
    v,q = TestFunctions(Z)
    if n_glen == 1.0:
        F = ( inner(Bn * D(u), D(v)) - p * div(v) - div(u) * q \
              - inner(f_body, v) ) * dx \
            - inner(outflow_sigma, v) * ds(bdryids['outflow'])
    else:
        Du2 = 0.5 * inner(D(u), D(u)) + (eps * Dtyp)**2.0
        rr = 1.0/n_glen - 1.0
        F = ( inner(Bn * Du2**(rr/2.0) * D(u), D(v)) - p * div(v) - div(u) * q - inner(f_body, v) ) * dx \
            - inner(outflow_sigma, v) * ds(bdryids['outflow'])

    # Dirichlet boundary conditions
    noslip = Constant((0.0, 0.0))
    inflow_u = getinflow(mesh,hsurfin,Hin,n_glen,alpha)
    bcs = [ DirichletBC(Z.sub(0), noslip, bdryids['base']),
            DirichletBC(Z.sub(0), inflow_u, bdryids['inflow']) ]

    # solve
    solve(F == 0, up, bcs=bcs,
          options_prefix='s',
          solver_parameters={"snes_converged_reason": True,
                             #"ksp_converged_reason": True,
                             #"ksp_monitor": True,
                             "ksp_type": "fgmres",  # or "gmres" or "minres"
                             "pc_type": "fieldsplit",
                             "pc_fieldsplit_type": "schur",
                             "pc_fieldsplit_schur_factorization_type": "full",  # or "diag"
                             "fieldsplit_0_ksp_type": "preonly",
                             "fieldsplit_0_pc_type": "lu",  # uses mumps in parallel
                             #"fieldsplit_0_ksp_converged_reason": True,
                             #"fieldsplit_1_ksp_converged_reason": True,
                             "fieldsplit_1_ksp_rtol": 1.0e-3,
                             "fieldsplit_1_ksp_type": "gmres",
                             "fieldsplit_1_pc_type": "none"})
    # ALSO CONSIDER:
    #    "ksp_type": "minres", "pc_type": "jacobi",
    #    "mat_type": "aij", "ksp_type": "preonly", "pc_type": "svd",  # fully-direct solver
    #    "mat_type": "aij", "ksp_view_mat": ":foo.m:ascii_matlab"
    # in parallel:  -s_fieldsplit_0_ksp_type gmres -s_fieldsplit_0_pc_type asm -s_fieldsplit_0_sub_pc_type ilu
    return up

