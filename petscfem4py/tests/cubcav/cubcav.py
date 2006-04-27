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

USE_SCHUR = PETSc.Options.HasName('use_schur')
USE_MATIS = PETSc.Options.HasName('use_matis')

petscopts = PETSc.Options()

if USE_SCHUR:
    USE_MATIS = True
    petscopts['ksp_type'] = 'preonly'
    petscopts['pc_type']  = 'schur'
    petscopts['pc_schur_ksp_type'] = 'gmres'
    petscopts['pc_schur_ksp_max_it'] = '250'
    petscopts['pc_schur_ksp_gmres_restart'] = '250'
    petscopts['pc_schur_pc_type']  = 'jacobi'
else:
    petscopts['ksp_type'] = 'gmres'
    if USE_MATIS:
        petscopts['pc_type'] = 'jacobi'
    else:
        petscopts['pc_type'] = 'bjacobi'


# Structured Grid
# ---------------

M,N,O = 4, 4, 4
M,N,O = 5, 5, 5
M,N,O = 10, 10, 10
M,N,O = 20, 20, 20
M,N,O = 30, 30, 30
M,N,O = 40, 40, 40
M,N,O = 50, 50, 50
#M,N,O = 60, 60, 60
#M,N,O = 70, 70, 70

nodes = numpy.arange(M*N*O)
nodes.shape = M,N,O

# Coordinates
# -----------
X, Y, Z = numpy.mgrid[0:1:M*1j, 0:1:N*1j, 0:1:O*1j]
xnod = numpy.array([X.ravel(), Y.ravel(), Z.ravel()])
xnod = xnod.transpose()


# Connectivity (Hexahedra)
# ------------
n0, n1, n2 = nodes.shape
e0, e1, e2 = (n-1 for n in (n0, n1, n2))
icone = numpy.zeros((e0, e1, e2, 8))
# front face
icone[...,0] = nodes[ 0:n0-1 , 0:n1-1, 0:n2-1 ]
icone[...,1] = nodes[ 0:n0-1 , 1:n1  , 0:n2-1 ]
icone[...,2] = nodes[ 0:n0-1 , 1:n1  , 1:n2   ]
icone[...,3] = nodes[ 0:n0-1 , 0:n1-1, 1:n2   ]
# back face
icone[...,4] = nodes[ 1:n0   , 0:n1-1, 0:n2-1 ]
icone[...,5] = nodes[ 1:n0   , 1:n1  , 0:n2-1 ]
icone[...,6] = nodes[ 1:n0   , 1:n1  , 1:n2   ]
icone[...,7] = nodes[ 1:n0   , 0:n1-1, 1:n2   ]
# final mesh
icone.shape = (e0*e1*e2, 8)


nnod = M*N*O
ndim = 3
ndof = ndim+1

nodeset = Nodeset(ndim, xnod)

elemset = Elemset('nsi_tet_les_fm2', icone)
elemset.setName('cubcav')
opts ={'geometry'   : 'cartesian%dd' % ndim,
       'ndim'       : str(ndim),
       'npg'        : str(8),
       'viscosity'  : str(1),
       'block_uploading': str(1),
       }
elemset.setOptions(opts)

# Boundary Conditions
# -------------------

dofset = Dofset(nnod, ndof)

# pressure
dofset.addFixations(0, 3, 0) # p

# top
nods = nodes[...,-1]
nods = nods.flatten()
_0, _1, _2 = [numpy.array(i).repeat(nods.size) for i in range(3)]
dofset.addFixations(nods, _0,  _1) # v_x = 1
dofset.addFixations(nods, _1,  _0) # v_y = 0
dofset.addFixations(nods, _2,  _0) # v_z = 0

# walls
walls = (nodes[:,:,0],
         nodes[0,:,:], nodes[-1,:,:],
         nodes[:,0,:], nodes[:,-1,:])
for nods in walls:
    nods = nods.flatten()
    _0, _1, _2 = [numpy.array(i).repeat(nods.size) for i in range(3)]
    dofset.addFixations(nods, _0,  _0) # v_x = 0
    dofset.addFixations(nods, _1,  _0) # v_y = 0
    dofset.addFixations(nods, _2,  _0) # v_z = 0


class NSI(NavierStokes):

    def __init__(self, nodeset, elemsetlist, dofset):

        fem = NavierStokes(nodeset, elemsetlist, dofset)

        COMM  = fem.getComm()
        sizes = fem.getDofSizes()
        
        J = PETSc.Mat()
        J.create(comm=COMM.comm)
        J.setSizes(sizes)
        if USE_MATIS:
            ldofs = fem.getLocalDofs()
            lgmap = PETSc.LGMapping(ldofs, comm=COMM.comm)
            J.setType(J.Type.IS)
            J.setLGMapping(lgmap)
            J.setPreallocation(27*4)
        else:
            J.setType(J.Type.MPIAIJ)
            J.setPreallocation([27*4, 9*4])
            
        x, r = J.getVecs()

        snes = PETSc.SNES(comm=J.comm)
        snes.setFunction(self.snes_res, r)
        snes.setJacobian(self.snes_jac, J)
        snes.setFromOptions()

        self.ass_time = 0.0
        self.sol_time = 0.0

        self.fem = fem
        self.J = J
        self.x = x
        self.r = r
        self.snes = snes
        self.t= 0

    def clear(self,*targs,**kargs):
        del self.fem
        del self.J
        del self.x
        del self.r
        del self.snes

    def snes_res(self, SNES, x, r):
        #PETSc.Print('computing res & jac...\n', SNES.comm)
        J,_,_ = SNES.getJacobian()

        r.zeroEntries()
        J.zeroEntries()
        
        t1 = wtime()
        self.fem.assemble(x, self.t, r, J)
        t2 = wtime()
        
        self.ass_time += t2-t1
        
    def snes_jac(self, SNES, x, J, P):
        assert (self.J == J)
        assert (self.J == P)
        if not J.isAssembled():
            J.assemble()
        return PETSc.Mat.Structure.SAME_NONZERO_PATTERN

    def solve(self, x, t=0.0):

        self.ass_time = 0.0
        self.sol_time = 0.0
        
        self.t  = t
        t1 = wtime()
        self.snes.solve(x)
        t2 = wtime()

        self.sol_time += t2-t1


nsi = NSI(nodeset, [elemset], dofset)
import atexit
atexit.register(nsi.clear)


if 0:
    # assemble performance
    # --------------------
    loops = 10
    snes = nsi.snes
    x = nsi.x
    r = x.duplicate()
    nsi.snes_res(snes,x,r)
    nsi.ass_time = 0.0
    for i in xrange(10):
        PETSc.Print('assemble loop: %d\n' % i)
        x.set(0)
        nsi.snes_res(snes,x,r)
    ass_time = PETSc.GlobalSum(nsi.ass_time)/SIZE
    PETSc.Print('Assemble time:  %f (%d loops)\n'% (ass_time, loops))
    PETSc.Print('\n')

if 1:
    PETSc.Print('Start solving...\n')
    x = nsi.x
    nsi.solve(x)
    ass_time = nsi.ass_time
    ass_time = PETSc.GlobalSum(ass_time)/SIZE
    sol_time = nsi.sol_time
    sol_time = PETSc.GlobalSum(sol_time)/SIZE
    snes_its = nsi.snes.iternum
    ksp_its  = nsi.snes.linear_its
    PETSc.Print('Assemble time:   %f\n'  % ass_time)
    PETSc.Print('Solution time:   %f\n'  % sol_time)
    PETSc.Print('SNES Iterations: %d\n'  % snes_its)
    PETSc.Print('KSP Iterations:  %d\n'  % ksp_its)


## dofmap = nsi.getDofMap()
## stt = numpy.array(state)
## sol = numpy.empty(nsi.getSize(), dtype=float)
## dofmap.apply(stt, sol, time[-1])
## sol.shape = (sol.size/ndof, ndof)


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
    
