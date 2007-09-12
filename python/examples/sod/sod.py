if __name__ == '__main__':
    import sys, petsc4py
    if 'petsc4py.PETSc' not in sys.modules:
        petsc4py.init(sys.argv)
    del sys, petsc4py

import numpy
from petsc4py import PETSc
from pf4py.kernel import *

opts = PETSc.Options()

D   = opts.getInt('D', 2)
N   = opts.getInt('N', 200)
L   = opts.getScalar('L', 10.0)

x0 = L/2.

xnod         = numpy.zeros((2*(N+1),2))
xnod[:,0]    = range(N+1)*2
xnod[N+1:,1] = 1.
xnod *= L/N

icone = numpy.array((range(N), range(1, N+1), range(N+2, 2*(N+1)), range(N+1, 2*N+1)))
icone = numpy.transpose(icone)

nodedata = DTableS(xnod)
elemconn = DTableI(icone+1)

nnod, ndim = nodedata.getShape()
ndof = ndim+2

domain = Domain('AD', ndim, nnod, ndof)

domain.setNodedata(nodedata)


R = opts.getReal('R', 287)

#elemset = Elemset('nsi_tet_les_full', elemconn)
elemset = Elemset('gasflow', elemconn)
elemset.setOptions({
    #'chunk_size'     : str(100000),
    #'block_uploading': str(1),
    
    'ndim'       : str(ndim),
    'geometry'   : 'cartesian%dd' % ndim,
    'npg'        : str(2**ndim),

    #'LES'    : str(0),
    #'C_smag' : str(0.18),

    'Rgas'   : str(287.),
    'ga'     : str(1.4),
    'visco'  : str(1e-15),
    'cond'   : str(1.0045e-12),

    'weak_form' : str(0),

    })

domain.addElemset(elemset)

# boundary conditions
# + velocity
nods_bot = numpy.arange(N+1)
nods_top = numpy.arange(N+1, 2*(N+1))

domain.setFixation(nods_bot, 2, 0.) # vy = 0
domain.setFixation(nods_top, 2, 0.)

#for i in range(ndof):
#    domain.setPeriodic(nods_bot, nods_top, i)

domain.setUp()

from pf4py.solvers import TransientSolver as Solver


solver = Solver(domain)

ts   = solver.ts
snes = ts.snes
ksp  = snes.ksp
pc   = ksp.pc

snes.rtol  = 1e-6
snes.atol  = 1e-12

T0 = opts.getReal('T0', 0.00)
dT = opts.getReal('dT', 0.0001)
nT = opts.getInt ('nT', 10)
T = T0 + nT*dT*1.001

ts.setTime(T0)
ts.setTimeStep(dT)
ts.setDuration(T, nT)

ts.setFromOptions()

#J = solver.jacobian
#J.convert('aijmumps', J)

#snes.setUseEW(True, version=3)

#ts.setUseFDColoring(True)


rhoL = opts.getReal('rhoL', 1.0)
rhoR = opts.getReal('rhoR', 0.125)
pL   = opts.getReal('pL',   1e5)
pR   = opts.getReal('pR',   1e4)

uini = numpy.zeros([nnod, ndof])
uini.shape = (nnod, -1)
uini[:, 0]      = numpy.where(xnod[:,0] <= x0, rhoL, rhoR)
uini[:, 1:ndof] = 0.
uini[:, ndof-1] = numpy.where(xnod[:,0] <= x0, pL,   pR)


try:
    sol = solver.run(uini)
except RuntimeError:
    PETSc.Error.view()
    raise

sol = sol.copy()

x   = xnod[0:N,0].ravel()
rho = sol[0:N,0].ravel()
u   = sol[0:N,1].ravel()
p   = sol[0:N,3].ravel()


if PETSc.COMM_WORLD.getSize() > 1: raise SystemExit

try:
    from matplotlib import pylab
except ImportError:
    print "'matplotlib' not available"
    raise SystemExit

pylab.figure()
pylab.plot(x,rho)

pylab.figure()
pylab.plot(x,u)

pylab.figure()
pylab.plot(x,p)

pylab.show()
