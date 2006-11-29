if __name__ == '__main__':
    import sys, petsc
    petsc.Init(sys.argv)
    del sys, petsc

import numpy
import petsc.PETSc as PETSc
from petscfem4py import *
from time import time as wtime

COMM = PETSc.COMM_WORLD
SIZE = COMM.size
RANK = COMM.rank

opts = PETSc.Options()

DUMP      = opts.hasName('dump')
USE_SCHUR = opts.hasName('schur')
USE_MATIS = opts.hasName('matis')

if USE_SCHUR and  SIZE > 1:
    USE_MATIS = True

PETSc.Sleep(2)

opts = PETSc.Options()

#opts['snes_max_it'] = '3'

if USE_SCHUR:
    opts['ksp_type'] = 'preonly'
    opts['pc_type']  = 'schur'
    #opts['pc_schur_ksp_type']         = 'gmres'
    #opts['pc_schur_ksp_left_pc']      = 'on'
    opts['pc_schur_ksp_right_pc']      = 'on'
    opts['pc_schur_ksp_max_it']        = '500'
    opts['pc_schur_ksp_gmres_restart'] = '500'
    opts['pc_schur_ksp_gmres_modifiedgramschmidt']= 'on'
else:
    #opts['ksp_type']         = 'gmres'
    #opts['ksp_left_pc']       = 'on'
    opts['ksp_right_pc']      = 'on'
    opts['ksp_max_it']        = '250'
    opts['ksp_gmres_restart'] = '250'
    opts['ksp_gmres_modifiedgramschmidt']= 'on'
    if USE_MATIS:
        opts['pc_type']     = 'jacobi'


grid_size = opts.getInt('grid', 10)
M,N = grid_size, grid_size

X, Y = numpy.mgrid[0:1:M*1j, 0:1:N*1j]
X, Y = X.transpose(), Y.transpose()
xnod = numpy.array([X.ravel(), Y.ravel()])
xnod = xnod.transpose()

if DUMP and RANK==0:
    fxnod = open('xnod.dat', 'w')
    for xyz in xnod:
        xyz.tofile(fxnod, sep=' ')
        fxnod.write('\n')
    fxnod.close()

nodes = numpy.arange(M*N).reshape(M,N)
nodes = numpy.flipud(nodes)
icone = numpy.zeros((M-1, N-1, 4))
icone[:,:,0] = nodes[ 1:M   , 0:N-1 ]
icone[:,:,1] = nodes[ 1:M   , 1:N   ]
icone[:,:,2] = nodes[ 0:M-1 , 1:N   ]
icone[:,:,3] = nodes[ 0:M-1 , 0:N-1 ]
icone.shape = (icone.size/4, 4)

if DUMP and RANK==0:
    ficone = open('icone.dat', 'w')
    for elem in (icone+1):
        elem.tofile(ficone, sep=' ')
        ficone.write('\n')
    ficone.close()

nnod = M*N
ndim = 2
ndof = 3

nodeset = Nodeset(ndim, xnod)

elemset = Elemset('nsi_tet_les_fm2', icone)
elemset.setName('fluid')
opt ={'chunk_size'     : str(40000),
      'block_uploading': str(1),
      'geometry'   : 'cartesian2d',
      'ndim'       : str(ndim),
      'npg'        : str(4),
      'rho'        : str(1),
      'viscosity'  : str(0.01),
      'viscosity'  : str(0.005),
      #'viscosity'  : str(0.001),
      }
elemset.setOptions(opt)


dofset = Dofset(nnod, ndof)

if DUMP and RANK==0:
    ffixa = open('fixa.dat', 'w')

if DUMP and RANK==0:
    ffixa.write('1 3 0.0\n')

# left, bottom, right
nods = numpy.r_[nodes[:,0], nodes[-1,:], nodes[:,-1]]
zero = numpy.zeros(nods.size)
one  = numpy.ones(nods.size)
dofset.addFixations(nods, zero, zero)
dofset.addFixations(nods, one,  zero)
dofset.addFixations(0, 2, 0)

if DUMP and RANK==0:
    for n in (nods+1):
        ffixa.write('%d %d %f\n' % (n, 1, 0.0))
        ffixa.write('%d %d %f\n' % (n, 2, 0.0))

# top
nods = nodes[0,:]
zero = numpy.zeros(nods.size)
one  = numpy.ones(nods.size)
two  = zero+2
dofset.addFixations(nods, zero, one)
dofset.addFixations(nods, one,  zero)

if DUMP and RANK==0:
    for n in (nods+1):
        ffixa.write('%d %d %f\n' % (n, 1, 1.0))
        ffixa.write('%d %d %f\n' % (n, 2, 0.0))

if DUMP and RANK==0:
    ffixa.close()

if DUMP:
    raise SystemExit


class System:

    def __init__(self, problem):

        self.problem = problem
        comm  = problem.getComm()
        ndofs = problem.getDofSizes()

        if USE_MATIS:
            ldofs = problem.getLocalDofs()
            lgmap = PETSc.LGMapping(ldofs, comm=comm.comm)
            J = PETSc.MatIS(ndofs, lgmap,  comm=comm.comm)
            J.getLocalMat().setType(J.Type.SEQAIJ)
            J.getLocalMat().setPreallocation(9*4)
        else:
            J = PETSc.Mat()
            J.create(comm=comm.comm)
            J.setSizes(ndofs)
            J.setType(J.Type.AIJ)
            J.setPreallocation([9*3, 3*3])

        x, r = J.getVecs()
        self.J, self.x0, self.r = J, x, r
        
        SNES = PETSc.SNES(comm=J.comm)
        SNES.setFunction(self.fun, r)
        SNES.setJacobian(self.jac, J)
        SNES.setFromOptions()
        self.SNES = SNES

    def fun(self, SNES, x, r):
        x0 = self.x0
        t0 = self.t0
        t  = self.t
        J,_,_ = SNES.getJacobian()
        r.zeroEntries()
        J.zeroEntries()
        t1 = wtime()
        #self.problem.assemble(x, t, r, J)
        #self.problem.assemble(x0, t0, x, t, r, J, self.alpha)
        self.problem.assemble(x0, t0, x, t, r, J, self.alpha, self.steady)
        t2 = wtime()
        self.assemble_time += t2-t1
        
    def jac(self, SNES, x, J, P):
        assert (self.J == J)
        assert (self.J == P)
        x0 = self.x0
        t0 = self.t0
        t  = self.t
        r,_ = SNES.getFunction()
        #r.zeroEntries()
        #J.zeroEntries()
        #self.problem.assemble(x0, t0, x, t, r, J, self.alpha, self.steady)

        return PETSc.Mat.Structure.SAME

    def step(self, x, t0, t, alpha=1.0, steady=False):
        self.t0 = t0
        self.t  = t
        self.alpha  = alpha
        self.steady = steady
        self.assemble_time = 0.0
        self.SNES.solve(x)

#alpha = 0.5
alpha  = 1.0
steady = True

problem = NavierStokes(nodeset, [elemset], dofset)

sys = System(problem)
state = sys.x0.duplicate()

#time = numpy.linspace(0, 100, 10)
#time = numpy.linspace(0, 100, 2)
time = numpy.linspace(0, 0.4, 5)


for i in xrange(1,len(time)):
    PETSc.Print('step: %d, time: %f\n' % (i, time[i]) )
    state.copy(sys.x0)
    sys.step(state, time[i-1], time[i], alpha, steady)
    x = sys.x0
    x.axpy(-1, state)
    PETSc.Print('du: %f\n' % x.norm())
    PETSc.Print('\n')

PETSc.Print('assemble time: %f\n'% sys.assemble_time, COMM)

#sys.SNES.view()

#PETSc.Print('\n')
#PETSc.SyncPrint("[%d] dof's: %d\n" % (RANK, sys.J.sizes[0][0]))
#PETSc.SyncFlush()

## dofmap = nsi.getDofMap()
## stt = numpy.array(state)
## sol = numpy.empty(nsi.getSize(), dtype=float)
## dofmap.apply(stt, sol, time[-1])
## sol.shape = (sol.size/3, 3)


## out = file('xnod.out', 'w')
## for xy in xnod:
##     out.write('%f %f\n' % tuple(xy))
## out.close()

## out = file('icone.out', 'w')
## for ijkl in icone:
##     out.write('%d %d %d %d\n' % tuple(ijkl))
## out.close()

## out = file('state.out', 'w')
## for ijk in sol:
##     out.write('%f %f %f\n' % tuple(ijk))
## out.close()

## if 0:
##     dofmap = nsi.getDofMap()

##     epart = elemset.getPart()
##     ldofs = nsi.getLocalDofs()

##     PETSc.Print('nodes:\n', COMM)
##     PETSc.Print('%s\n' % nodes, COMM)

##     PETSc.Print('icone:\n', COMM)
##     PETSc.Print('%s\n' % icone, COMM)

##     PETSc.Print('partitioning:\n', COMM)
##     PETSc.SyncPrint('[%d]: %s\n' % (RANK, list(epart)), COMM)
##     PETSc.SyncFlush(COMM)
##     PETSc.Print('\n', COMM)

##     PETSc.Print('local dofs:\n', COMM)
##     PETSc.SyncPrint('[%d]: %d %s\n' % (RANK, ldofs.size, list(ldofs)), COMM)
##     PETSc.SyncFlush(COMM)
##     PETSc.Print('\n', COMM)

## J = sys.J
## J.assemble()
## J.zeroEntries()
## J.view()

## ldofs = nsi.getLocalDofs()
## lgmap = PETSc.LGMapping(ldofs)
## gdofs = lgmap.applyInverse(ldofs)

## lnods = ldofs[::ndof]/ndof
## PETSc.Print('local nodes:\n', COMM)
## PETSc.SyncPrint('[%d]: %d %s\n' % (RANK, lnods.size, list(lnods)), COMM)
## PETSc.SyncFlush(COMM)
## PETSc.Print('\n', COMM)
