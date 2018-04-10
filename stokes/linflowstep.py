#!/usr/bin/env python3
# (C) 2018 Ed Bueler

# Solve Newtonian fluid glacier bedrock step Stokes problem.
# See mccarthy/projects/2018/flowstep/README.md.

# Usage with default bedstep:
#     $ ./genstepmesh.py glacier.geo      # create domain
#     $ gmsh -2 glacier.geo               # mesh domain
#     $ source ~/firedrake/bin/activate
#     (firedrake) $ ./linflowstep.py      # solve Stokes problem
#     (firedrake) $ paraview glacier.pvd  # visualize

# Note that the genstepmesh.py script allows uniform refinement
# (-refine X) and refinement at interior corner (-refine_factor X).

# Usage with zero bedstep for slab-on-slope:
#     $ ./genstepmesh.py -bs 0.0 slab.geo
#     $ gmsh -2 slab.geo
#     $ source ~/firedrake/bin/activate
#     (firedrake) $ ./linflowstep.py -bs 0.0 -f slab

import argparse

# process options
defaultmix = 'P2P1'
mixchoices = ['P2P1','P3P2','P2P0','CRP0','P1P0']
defaultf = 'glacier'
parser = argparse.ArgumentParser(description='Solve Newtonian-fluid glacier bedrock step Stokes problem.')
parser.add_argument('-elements', metavar='X', default=defaultmix, choices=mixchoices,
                    help='stable mixed finite elements from: %s (default=%s)' % (','.join(mixchoices),defaultmix) )
parser.add_argument('-f', metavar='ROOT', default=defaultf,
                    help='input/output file name root (default=%s)' % defaultf)
parser.add_argument('-bs', type=float, default=120.0, metavar='BS',
                    help='height of bed step (m; default = 120.0)')
args, unknown = parser.parse_known_args()
inname = args.f + '.msh'
outname = args.f + '.pvd'
bs = args.bs

from firedrake import *

# glacier physical constants
g = 9.81                   # m s-2
rho = 910.0                # kg m-3
A_ice = 3.1689e-24         # Pa-3 s-1; EISMINT I value of ice softness
B_ice = A_ice**(-1.0/3.0)  # Pa s0.33333;  ice hardness
secpera = 31556926.0       # seconds per year
D_typical = 10.0 / secpera
nu_e = B_ice * D_typical**(-2.0/3.0)  # effective viscosity FIXME for Glen law
print('effective viscosity %.2e Pa s' % nu_e)

# input mesh and define geometry
print('reading mesh from %s ...' % inname)
mesh = Mesh(inname)
print('mesh has %d vertices and %d elements' % (mesh.num_vertices(),mesh.num_cells()))
top_id = 32
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
u,p = TrialFunctions(Z)
v,q = TestFunctions(Z)

# define weak form
a = ( nu_e * inner(grad(u), grad(v)) - p * div(v) - div(u) * q ) * dx

# alternative from Mitchell thesis (2012) seems to cause pressure
# concentration at bottom of outflow:
#a = ( 0.5 * nu_e * inner(grad(v)+grad(v).T, grad(u)+grad(u).T)
#      - p * div(v) - div(u) * q ) * dx

# define body force
f = Constant((g * rho * sin(alpha), - g * rho * cos(alpha)))

# right side includes hydrostatic normal force on outflow
x,z = SpatialCoordinate(mesh)
outflow_sigma = as_vector([- rho * g * cos(alpha) * (H - z), 0.0])
L = inner(f, v) * dx + inner(outflow_sigma, v) * ds(outflow_id)

# Dirichlet boundary conditions, defined on the velocity space, apply
# on base and inflow
noslip = Constant((0.0, 0.0))
# assume Newtonian slab-on-slope inflow (uses nu_e)
C = rho * g * sin(alpha) / nu_e
inflow_u = as_vector([C * ((H + bs) * (z - bs) + (bs*bs - z*z)/2.0), 0.0])
bcs = [ DirichletBC(Z.sub(0), noslip, (top_id,)),
        DirichletBC(Z.sub(0), inflow_u, (inflow_id,)) ]

# solve
print('solving linear variational problem with %s elements ...' % args.elements)
up = Function(Z)   # put solution here
solve(a == L, up, bcs=bcs,
      options_prefix='lfs',
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

# note default solver already has -snes_type ksponly for this linear problem
# ALSO:
#    "ksp_type": "minres", "pc_type": "jacobi",
#    "mat_type": "aij", "ksp_type": "preonly", "pc_type": "svd",  # fully-direct solver
#    "mat_type": "aij", "ksp_view_mat": ":foo.m:ascii_matlab"

u,p = up.split()
u.rename('velocity')
p.rename('pressure')

one = Function(W)
one.interpolate(0.0*one + 1.0)
area = assemble(dot(one,one) * dx)
print('domain area = %.2e m2' % area)
pav = assemble(sqrt(dot(p, p)) * dx) / area
print('average pressure = %.2f Pa' % pav)
velmagav = assemble(sqrt(dot(u, u)) * dx) / area
print('average velocity magnitude = %.2f m a-1' % (secpera * velmagav))

# write the solution to a file for visualisation with paraview
print('writing (velocity,pressure) to %s ...' % outname)
File(outname).write(u,p)

