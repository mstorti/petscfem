if __name__ == '__main__':
    import sys, petsc
    petsc.Init(sys.argv)
    del sys, petsc

import numpy
import petsc.PETSc as PETSc
from petscfem4py import *
from time import time as wtime

USE_SCHUR = PETSc.Options.HasName('use_schur')
USE_MATIS = PETSc.Options.HasName('use_matis')
if USE_SCHUR:
    USE_MATIS = True
    
COMM = PETSc.COMM_WORLD
SIZE = COMM.size
RANK = COMM.rank

petscopts = PETSc.Options()

if USE_MATIS and USE_SCHUR:
    petscopts['ksp_type']     = 'preonly'
    petscopts['ksp_max_it']   = 1

    petscopts['pc_type']      = 'schur'
    petscopts['schur_global_ksp_type']    = 'gmres'
    petscopts['schur_global_pc_type']     = 'jacobi'
    #petscopts['schur_global_ksp_monitor'] = 'stdout'
else:
    petscopts['ksp_type']    = 'gmres'
    petscopts['pc_type']     = 'jacobi'
    #petscopts['ksp_monitor'] = 'stdout'


petscopts['snes_monitor'] = 'stdout'
petscopts['snes_rtol'] = '1e-5'




M, N = 15, 15

X, Y = numpy.mgrid[0:1:M*1j, 0:1:N*1j]
X, Y = X.transpose(), Y.transpose()
xnod = numpy.array([X.ravel(), Y.ravel()])
xnod = xnod.transpose()

nodes = numpy.arange(M*N).reshape(M,N)
nodes = numpy.flipud(nodes)
icone = numpy.zeros((M-1, N-1, 4))
icone[:,:,0] = nodes[ 1:M   , 0:N-1 ]
icone[:,:,1] = nodes[ 1:M   , 1:N   ]
icone[:,:,2] = nodes[ 0:M-1 , 1:N   ]
icone[:,:,3] = nodes[ 0:M-1 , 0:N-1 ]
icone.shape = (icone.size/4, 4)


nnod = M*N
ndim = 2
ndof = 3

nodeset = Nodeset(ndim, xnod)

elemset = Elemset('nsi_tet_les_fm2', icone)
elemset.setName('malla')
opts ={'geometry'   : 'cartesian2d',
       'ndim'       : str(ndim),
       'npg'        : str(4),
       #'viscosity'  : str(0.001),
       'viscosity'  : str(1),
       }
elemset.setOptions(opts)


dofset = Dofset(nnod, ndof)

# left, bottom, right
nods = numpy.r_[nodes[:,0], nodes[-1,:], nodes[:,-1]]
zero = numpy.zeros(nods.size)
one  = numpy.ones(nods.size)
dofset.addFixations(nods, zero, zero)
dofset.addFixations(nods, one,  zero)
dofset.addFixations(0, 2, 0)

# top
nods = nodes[0,:]
zero = numpy.zeros(nods.size)
one  = numpy.ones(nods.size)
two  = zero+2
dofset.addFixations(nods, zero, one)
dofset.addFixations(nods, one,  zero)




nsi = NavierStokes(nodeset, [elemset], dofset)

class System:

    def __init__(self, problem):
        self.problem = problem

        comm  = problem.getComm()
        ndofs = problem.getDofSizes()


        if USE_MATIS:
            ldofs = nsi.getLocalDofs()
            lgmap = PETSc.LGMapping(ldofs, comm=comm.comm)
            J = PETSc.MatIS(ndofs, lgmap,  comm=comm.comm)
            J.getLocalMat().setType(J.Type.SEQAIJ)
            J.getLocalMat().setPreallocation(9*4)
        else:
            J = PETSc.Mat()
            J.create(comm=comm.comm)
            J.setSizes(ndofs)
            J.setType(J.Type.AIJ)
            J.setPreallocation(9*4)

        x, r = J.getVecs()
        self.J, self.x, self.r = J, x, r
        
        SNES = PETSc.SNES(comm=J.comm)
        SNES.setFunction(self.fun, r)
        SNES.setJacobian(self.jac, J)
        SNES.setFromOptions()
        self.SNES = SNES

    def fun(self, SNES, x, r):
        xi = self.x
        ti = self.ti
        t  = self.t
        J,_,_ = SNES.getJacobian()
        r.zeroEntries()
        J.zeroEntries()
        #self.problem.assemble(xi, ti, x, t, r, J, self.alpha)
        t1 = wtime()
        self.problem.assemble(x, t, r, J)
        t2 = wtime()
        self.assemble_time += t2-t1
        
    def jac(self, SNES, x, J, P):
        assert (self.J == J)
        assert (self.J == P)
        #J.assemble()
        return PETSc.Mat.Structure.SAME

    def step(self, x, ti, t, alpha=1.0):
        x.copy(self.x)
        self.ti = ti
        self.t  = t
        self.alpha = alpha
        self.assemble_time = 0.0
        self.SNES.solve(x)


sys = System(nsi)
state = sys.x.duplicate()
#time = numpy.linspace(0, 100, 10)
time = numpy.linspace(0, 100, 2)

#time = []

#alpha = 0.5
alpha = 1
for i in xrange(len(time)-1):
    PETSc.Print('step: %d, time: %f\n' % (i, time[i+1]) )
    sys.step(state, time[i], time[i+1], alpha)
    sys.x.axpy(-1,state)
    PETSc.Print('du: %f\n' % sys.x.norm())
    PETSc.Print('\n')

PETSc.Print('assemble time: %f\n'% sys.assemble_time, COMM)

#sys.SNES.view()

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

if 0:
    dofmap = nsi.getDofMap()

    epart = elemset.getPart()
    ldofs = nsi.getLocalDofs()

    PETSc.Print('nodes:\n', COMM)
    PETSc.Print('%s\n' % nodes, COMM)

    PETSc.Print('icone:\n', COMM)
    PETSc.Print('%s\n' % icone, COMM)

    PETSc.Print('partitioning:\n', COMM)
    PETSc.SyncPrint('[%d]: %s\n' % (RANK, list(epart)), COMM)
    PETSc.SyncFlush(COMM)
    PETSc.Print('\n', COMM)

    PETSc.Print('local dofs:\n', COMM)
    PETSc.SyncPrint('[%d]: %d %s\n' % (RANK, ldofs.size, list(ldofs)), COMM)
    PETSc.SyncFlush(COMM)
    PETSc.Print('\n', COMM)

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
