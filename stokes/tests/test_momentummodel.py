from firedrake import *
from momentummodel import secpera, dayspera, MomentumModel
from domain import bdryids

# see bottom page https://www.firedrakeproject.org/extruded-meshes.html
extrudebdryids = {'inflow'  : 1,
                  'outflow' : 2,
                  'top'     : 'top',
                  'base'    : 'bottom'}

def test_slab_extrude():
    L = 3000.0
    H0 = 400.0
    mx = 20
    mz = 5
    basemesh = IntervalMesh(mx, length_or_left=L)
    mesh = ExtrudedMesh(basemesh, layers=mz, layer_height=H0 / mz)
    mm = MomentumModel()
    up = mm.solve(mesh, extrudebdryids, extrudemode=True, package='Direct')
    umagav,umagmax,pav,pmax = mm.solutionstats(mesh)
    assert abs(secpera * umagav  - 724.5) < 1.0  # check speeds in m a-1
    assert abs(secpera * umagmax - 907.3) < 1.0
    assert abs(1.0e-5 * pav      - 17.8) < 0.1   # check pressures in bar
    assert abs(1.0e-5 * pmax     - 35.5) < 0.1
    #print('flow speed: av = %10.3f m a-1,  max = %10.3f m a-1' % (secpera*umagav,secpera*umagmax))
    #print('pressure:   av = %10.3f bar,    max = %10.3f bar' % (1.0e-5*pav,1.0e-5*pmax))
    #u, p = up.subfunctions
    #from firedrake.output import VTKFile
    #VTKFile('foo.pvd').write(u, p)

def test_bedstep_coarse():
    mesh = Mesh('coarse.msh')
    mm = MomentumModel()
    up = mm.solve(mesh, bdryids, package='Direct')
    umagav,umagmax,pav,pmax = mm.solutionstats(mesh)
    assert mesh.num_cells() == 180
    assert mesh.num_vertices() == 118
    assert abs(secpera * umagav  - 592.6) < 0.1  # check speeds in m a-1
    assert abs(secpera * umagmax - 906.1) < 0.1
    assert abs(1.0e-5 * pav      - 17.46) < 0.1   # check pressures in bar
    assert abs(1.0e-5 * pmax     - 40.06) < 0.1
    #print('mesh has %d elements (cells) and %d vertices' % (mesh.num_cells(),mesh.num_vertices()))
    #print('flow speed: av = %10.3f m a-1,  max = %10.3f m a-1' % (secpera*umagav,secpera*umagmax))
    #print('pressure:   av = %10.3f bar,    max = %10.3f bar' % (1.0e-5*pav,1.0e-5*pmax))
    #u, p = up.subfunctions
    #from firedrake.output import VTKFile
    #VTKFile('foo.pvd').write(u, p)

#test_slab_extrude()
#test_bedstep_coarse()