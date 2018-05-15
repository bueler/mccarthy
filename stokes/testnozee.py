import argparse
parser = argparse.ArgumentParser(\
    description='test idea about removing z-dependence of a scalar function')
parser.add_argument('inname', metavar='INNAME',
                    help='input file name ending with .msh')
args, unknown = parser.parse_known_args()

from firedrake import *

def nozee(mesh,X,top_id,bottom_id,Gtop):
    G = TrialFunction(X)
    v = TestFunction(X)
    khat = Constant((0.0,1.0))
    a = - G * inner(khat, grad(v)) * dx \
        + G * v * inner(khat,FacetNormal(mesh)) * ds(bottom_id)
    L = Constant(0.0) * v * dx
    bc = [ DirichletBC(X, Gtop, top_id) ]
    Gsoln = Function(X)
    solve(a == L, Gsoln, bcs=bc, options_prefix='nz',
          solver_parameters={"ksp_type": "preonly",
                             "pc_type": "svd"})
          # also -nz_pc_type lu -nz_pc_factor_nonzeros_along_diagonal
    return Gsoln

from gendomain import L,bdryids

mesh = Mesh(args.inname)
print('mesh in %s has %d elements (cells) and %d vertices' \
      % (args.inname,mesh.num_cells(),mesh.num_vertices()))

X = FunctionSpace(mesh, "CG", 2)
x,z = SpatialCoordinate(mesh)

Gtop = sin(6*pi*x/L)
G = nozee(mesh,X,bdryids['top'],bdryids['base'],Gtop)
G.rename('G')
File('nozee.pvd').write(G)

