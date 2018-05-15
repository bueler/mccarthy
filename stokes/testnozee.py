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
    return Gsoln

from gendomain import L,bdryids
mesh = Mesh('glacier.msh')
print('mesh in "glacier.msh" has %d elements (cells) and %d vertices' \
      % (mesh.num_cells(),mesh.num_vertices()))

X = FunctionSpace(mesh, "CG", 1)
x,z = SpatialCoordinate(mesh)

Gtop = sin(6*pi*x/L)
G = nozee(mesh,X,bdryids['top'],bdryids['base'],Gtop)
G.rename('G')
File('nozee.pvd').write(G)

