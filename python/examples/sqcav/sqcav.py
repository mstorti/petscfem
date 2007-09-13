if __name__ == '__main__':
    import sys, petsc4py
    petsc4py.init(sys.argv)
    del sys, petsc4py

import numpy
from petsc4py import PETSc
from pf4py.kernel import DTableS, DTableI, Domain, Elemset
from pf4py.solvers import TransientSolver as Solver
from pf4py.tools.boxmesh import BoxMesh

opts = PETSc.Options()


N = opts.getInt ('N',   32)
L = opts.getReal('L',   1.0)
h = L/(N-1)

U   = opts.getReal('U',   1.0)
rho = opts.getReal('rho', 1.0)
Re  = opts.getReal('Re',  100)

CFL = opts.getReal('CFL',  1.0)
dT  = opts.getReal('dT',   h*CFL/U)
T0  = opts.getReal('T0',   0.0)
nT  = opts.getInt ('nT',   5)

xnod, icone, wall = BoxMesh([N,N], bbox=[(0,L), (0,L)])

nnod, ndim = xnod.shape
ndof = ndim + 1 
domain = Domain('NS', ndim, nnod, ndof, PETSc.COMM_WORLD)
domain.setNodedata(DTableS(xnod))
elemset = Elemset('nsi_tet_les_full')
elemset.setData(DTableI(icone+1))
elemset.setOptions({
    'rho'        : str(rho),
    'viscosity'  : str(rho*U*L/Re),
    'use_lsic'   : str(1), # XXX
    'use_tst'    : str(1), # XXX
    })
domain.addElemset(elemset)
# boundary conditions
# + velocity
for k in wall:
    if k == 'top': # lid
        domain.setFixation(wall[k], 0, U)   # vx = 1
        domain.setFixation(wall[k], 1, 0.0) # vy = 0
    else:
        domain.setFixation(wall[k], 0, 0.0) # vx = 0
        domain.setFixation(wall[k], 1, 0.0) # vy = 0
# + pressure
domain.setFixation(0, 2, 0.0)               # p = 0
# XXX
domain.setUp()


solver = Solver(domain, alpha=0.55)

ts   = solver.ts
snes = ts.snes
ksp  = snes.ksp
pc   = ksp.pc


ts.time      = T0
ts.time_step = dT
ts.max_steps = nT
ts.max_time  = T0 + 1.001 * nT*dT

snes.rtol = 1e-6
snes.atol = 1e-12

ksp.rtol  = 1e-4
ksp.atol  = 1e-12

ksp  = snes.ksp


## def save_solution(ts, step, time, u):
##     domain = ts.getAppCtx()
##     sol = PETSc.Vec()
##     domain.allocateSolution(sol)
##     domain.buildSolution(time, u, sol)
##     asol = numpy.asarray(sol)
##     asol.shape = (domain.getNNod(), domain.getNDof())
##     asol.tofile('state_%d.out' % step)
## ts.setMonitor(save_solution)

ts.setFromOptions()

u0 = None
soln = solver.run(u0)

x, y = numpy.mgrid[0:1:N*1j,0:1:N*1j]
vx = soln[:,0].reshape(N,N)
vy = soln[:,1].reshape(N,N)
p  = soln[:,2].reshape(N,N)
u  = numpy.sqrt(vx**2 + vy**2)


if domain.getComm().getSize() > 1: raise SystemExit

try:
    from matplotlib import pylab
except ImportError:
    print "'matplotlib.pylab' not available"
    raise SystemExit

pylab.figure()
pylab.quiver(x,y,vx,vy,scale=5)
pylab.axis('equal')

pylab.figure()
pylab.contourf(x,y,u)
pylab.axis('equal')

pylab.figure()
cf = pylab.contourf(x,y,p)
pylab.axis('equal')

pylab.show()
