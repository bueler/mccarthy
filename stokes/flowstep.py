#!/usr/bin/env python3
# (C) 2018 Ed Bueler

# FIXME  note other FIXMEs below; this version sort of works with n=1:
# ./flowstep.py -bs 0.0 -n_glen 1.0 -f slab -s_snes_max_it 1000 -s_snes_rtol 1.0e-4

# Solve glacier bedrock-step Glen-Stokes problem.
# See mccarthy/projects/2018/flowstep/README.md.

# Default usage:
#     $ ./genstepmesh.py glacier.geo      # create domain
#     $ gmsh -2 glacier.geo               # mesh domain
#     $ source ~/firedrake/bin/activate
#     (firedrake) $ ./flowstep.py         # solve Stokes problem
#     (firedrake) $ paraview glacier.pvd  # visualize

# Note that the genstepmesh.py script allows uniform refinement
# (-refine X) and refinement at interior corner (-refine_corner X).

# Usage with zero bedstep for slab-on-slope:
#     $ ./genstepmesh.py -bs 0.0 slab.geo
#     $ gmsh -2 slab.geo
#     $ source ~/firedrake/bin/activate
#     (firedrake) $ ./flowstep.py -bs 0.0 -f slab

import argparse

# process options
defaultmix = 'P2P1'
mixchoices = ['P2P1','P3P2','P2P0','CRP0','P1P0']
defaultf = 'glacier'
parser = argparse.ArgumentParser(description='Solve glacier bedrock-step Glen-Stokes problem.')
parser.add_argument('-elements', metavar='X', default=defaultmix, choices=mixchoices,
                    help='stable mixed finite elements from: %s (default=%s)' % (','.join(mixchoices),defaultmix) )
parser.add_argument('-f', metavar='ROOT', default=defaultf,
                    help='input/output file name root (default=%s)' % defaultf)
parser.add_argument('-bs', type=float, default=120.0, metavar='X',
                    help='height of bed step (m; default = 120.0)')
parser.add_argument('-initonly', action='store_true')
parser.add_argument('-n_glen', type=float, default=3.0, metavar='X',
                    help='Glen flow law exponent (default = 3.0)')
args, unknown = parser.parse_known_args()
inname = args.f + '.msh'
outname = args.f + '.pvd'
bs = args.bs
n_glen = args.n_glen

from firedrake import *

# glacier physical constants
secpera = 31556926.0       # seconds per year
g = 9.81                   # m s-2
rho = 910.0                # kg m-3
A_ice = 3.1689e-24         # Pa-3 s-1; EISMINT I value of ice softness; tied to n=3
B_ice = A_ice**(-1.0/3.0)  # Pa s(1/3);  ice hardness

# input mesh and define geometry
print('reading mesh from %s ...' % inname)
mesh = Mesh(inname)
print('mesh has %d vertices and %d elements' % (mesh.num_vertices(),mesh.num_cells()))
base_id = 32
inflow_id = 33
outflow_id = 34
H = 400.0
alpha = 0.1

# define mixed finite element space
mixFE = {'P2P1'  : (VectorFunctionSpace(mesh, "CG", 2), # Taylor-Hood
                    FunctionSpace(mesh, "CG", 1)),
         'P3P2'  : (VectorFunctionSpace(mesh, "CG", 3),
                    FunctionSpace(mesh, "CG", 2)),
         'P2P0'  : (VectorFunctionSpace(mesh, "CG", 2),
                    FunctionSpace(mesh, "DG", 0)),
         'CRP0'  : (VectorFunctionSpace(mesh, "CR", 1),
                    FunctionSpace(mesh, "DG", 0)),
         'P1P0'  : (VectorFunctionSpace(mesh, "CG", 1), # interesting but not
                    FunctionSpace(mesh, "DG", 0))       #   recommended
        }
(V,W) = mixFE[args.elements]
Z = V * W
v,q = TestFunctions(Z)

# define body force
f_body = Constant((g * rho * sin(alpha), - g * rho * cos(alpha)))

# Dirichlet boundary condition applied on base
noslip = Constant((0.0, 0.0))

# right side outflow: apply hydrostatic normal force; nonhomogeneous Neumann
x,z = SpatialCoordinate(mesh)
outflow_sigma = as_vector([- rho * g * cos(alpha) * (H - z), 0.0])

# linear problem is how we generate initial iterate
def initialiterate(Z):
    # assume Newtonian slab-on-slope inflow
    D_typical = 10.0 / secpera
    nulin = B_ice * D_typical**(-2.0/3.0)  # effective viscosity
    print('effective viscosity in initialization %.2e Pa s' % nulin)
    Clin = rho * g * sin(alpha) / nulin
    inflowlin = as_vector([Clin * ((H + bs) * (z - bs) + (bs*bs - z*z)/2.0), 0.0])
    ulin,plin = TrialFunctions(Z)
    bcslin = [ DirichletBC(Z.sub(0), noslip, (base_id,)),
               DirichletBC(Z.sub(0), inflowlin, (inflow_id,)) ]
    # Newtonian weak form
    alin = ( nulin * inner(grad(ulin), grad(v)) - plin * div(v) - div(ulin) * q ) * dx
    Llin = inner(f_body, v) * dx + inner(outflow_sigma, v) * ds(outflow_id)
    print('solving linear variational problem for initial iterate ...')
    uplin = Function(Z)   # put solution here
    solve(alin == Llin, uplin, bcs=bcslin,
          options_prefix='lin',
          solver_parameters={#"ksp_view": True,
                             "ksp_converged_reason": True,
                             "ksp_monitor": True,
                             "ksp_type": "fgmres",  # or "gmres" or "minres"
                             "pc_type": "fieldsplit",
                             "pc_fieldsplit_type": "schur",
                             "pc_fieldsplit_schur_factorization_type": "full",  # or "diag"
                             "fieldsplit_0_ksp_type": "preonly",
                             "fieldsplit_0_pc_type": "lu",
                             #"fieldsplit_0_ksp_converged_reason": True,
                             "fieldsplit_1_ksp_converged_reason": True,
                             "fieldsplit_1_ksp_rtol": 1.0e-3,
                             "fieldsplit_1_ksp_type": "gmres",
                             "fieldsplit_1_pc_type": "none"})
    return uplin

# put solution here
up = Function(Z)       # *not* TrialFunctions(Z)
u,p = split(up)        # up.split() not equivalent here?

# define the nonlinear weak form F(u;v)
D_typical = 10.0 / secpera
eps2 = 0.0001 * D_typical**2.0
Du = 0.5 * (grad(u) + grad(u).T)
normsqrDu = 0.5 * inner(Du, Du) + eps2
rr = 0.5 * (1.0/n_glen - 1.0)

#F = ( inner(B_ice * normsqrDu**rr * Du, grad(v)) - p * div(v) - div(u) * q - inner(f_body, v) ) * dx \
#    - inner(outflow_sigma, v) * ds(outflow_id)

# FIXME:  choose rr=0  and  Du = grad(u)
F = ( inner(B_ice * grad(u), grad(v)) - p * div(v) - div(u) * q - inner(f_body, v) ) * dx \
    - inner(outflow_sigma, v) * ds(outflow_id)

# inflow boundary condition
Hin = H + bs
if False:
    # assume Newtonian slab-on-slope inflow
    nulin = B_ice * D_typical**(-2.0/3.0)  # effective viscosity
    Clin = rho * g * sin(alpha) / nulin
    inflow_u = as_vector([Clin * (Hin * (z - bs) + (bs*bs - z*z)/2.0), 0.0])
else:
    # assume slab-on-slope inflow
    C = 2.0 * A_ice * (rho * g * sin(alpha))**n_glen / (n_glen + 1.0)
    inflow_u = as_vector([C * (H**(n_glen+1.0) - (Hin - z)**(n_glen+1.0)), 0.0])

bcs = [ DirichletBC(Z.sub(0), noslip, (base_id,)),
        DirichletBC(Z.sub(0), inflow_u, (inflow_id,)) ]

# solve
up0 = initialiterate(Z)
u0,p0 = up0.split()
u,p = up.split()
u.interpolate(u0)
p.interpolate(p0)

if not args.initonly:
    print('solving nonlinear variational problem ...')
    solve(F == 0, up, bcs=bcs,
          options_prefix='s',
          solver_parameters={"snes_mf": True,  # FIXME only for very coarse grids
                             "mat_type": "aij",
                             "ksp_type": "preonly",
                             "pc_type": "svd",
                             #"snes_rtol": 1.0e-14,
                             #"snes_stol": 0.0,
                             "snes_monitor": True,
                             "snes_max_funcs": 100000,
                             "snes_converged_reason": True,
                             "ksp_converged_reason": True,
                             "ksp_monitor": True})

# see linflowstep.py for more solver combinations

u,p = up.split()
u.rename('velocity')
p.rename('pressure')
# examine values directly:   print(p.vector().array())

one = Constant(1.0, domain=mesh)
area = assemble(dot(one,one) * dx)
print('domain area = %.2e m2' % area)
pav = assemble(sqrt(dot(p, p)) * dx) / area
print('average pressure = %.2f Pa' % pav)
velmagav = assemble(sqrt(dot(u, u)) * dx) / area
print('average velocity magnitude = %.2f m a-1' % (secpera * velmagav))

# write the solution to a file for visualisation with paraview
print('writing (velocity,pressure) to %s ...' % outname)
File(outname).write(u,p)

