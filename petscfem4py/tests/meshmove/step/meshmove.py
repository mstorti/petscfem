import numpy
import petsc.PETSc as PETSc
from petscfem4py import *

M, N = 11, 11

X, Y = numpy.mgrid[0:1:M*1j, 0:1:N*1j]
X, Y = X.transpose(), Y.transpose()
xnod = numpy.array([X.ravel(), Y.ravel()])
xnod = xnod.transpose()

nodes = numpy.arange(M*N).reshape(M,N)
nodes = numpy.flipud(nodes)
icone = numpy.zeros((M-1, N-1, 4))
icone[:,:,0] = nodes[ 1:M   , 0:N-1   ]
icone[:,:,1] = nodes[ 1:M   , 1:N   ]
icone[:,:,2] = nodes[ 0:M-1 , 1:N   ]
icone[:,:,3] = nodes[ 0:M-1 , 0:N-1   ]
icone.shape = (icone.size/4, 4)
icone = numpy.r_[ icone[:, (0,1,2,)], icone[:, (2,3,0)]]


nnod = M*N
ndim = 2
ndof = 2

nodeset = Nodeset(ndim, xnod)

elemset = Elemset('mesh_move', icone)
elemset.setName('malla')
opts ={'geometry'   : 'triangle',
       'ndim'       : str(ndim),
       'npg'        : str(1),
       'distor_exp' : str(-1),
       }
elemset.setOptions(opts)


dofset = Dofset(nnod, ndof)

# time dependet boundary conditions
class Ramp(Amplitude):
    def __call__(self, node, field, time):
        return time
class Tanh(Amplitude):
    def __call__(self, node, field, time):
        return numpy.tanh(3*time)
amp = Ramp()

# top
nods = nodes[0,:]
zero = numpy.zeros(nods.size)
one  = numpy.ones(nods.size)
vals = numpy.r_[1:0:nods.size*1j]
dofset.addFixations(nods, zero, vals, amp)
dofset.addFixations(nods, one,  zero)

# left
half = M//2

nods = nodes[:half,0]
zero = numpy.zeros(nods.size)
one  = numpy.ones(nods.size)
dofset.addFixations(nods, zero, one, amp)

nods = nodes[half:,0]
zero = numpy.zeros(nods.size)
vals = numpy.r_[1:0:nods.size*1j]
dofset.addFixations(nods, zero, vals, amp)

nods = nodes[:,0]
zero = numpy.zeros(nods.size)
one  = numpy.ones(nods.size)
dofset.addFixations(nods, one,  zero)

# bottom
nods = nodes[-1,:]
zero = numpy.zeros(nods.size)
one  = numpy.ones(nods.size)
for field in [zero, one]:
    dofset.addFixations(nods, field, zero)

# right
nods = nodes[:,-1]
zero = numpy.zeros(nods.size)
one  = numpy.ones(nods.size)
for field in [zero, one]:
    dofset.addFixations(nods, field, zero)


nsi = NavierStokes(nodeset, [elemset], dofset)

class System:

    def __init__(self, problem):
        self.problem = problem

        comm  = problem.getComm()
        ndofs = problem.getDofSizes()
        J  = PETSc.MatSeqAIJ(ndofs, nz=9*4, comm=comm.comm)
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
        self.problem.assemble(xi, ti, x, t, r, J, self.alpha)
        
    def jac(self, SNES, x, J, P):
        assert (self.J == J)
        assert (self.J == P)
        return PETSc.Mat.Structure.SAME

    def step(self, x, ti, t, alpha=1.0):
        x.copy(self.x)
        self.ti = ti
        self.t  = t
        self.alpha = alpha
        self.SNES.solve(x)


petscopts = PETSc.Options()
petscopts['ksp_type']     = 'preonly'
petscopts['pc_type']      = 'lu'
petscopts['snes_monitor'] = 'stdout'
#petscopts['ksp_monitor']  = 'stdout'
#petscopts['snes_view']    = 'stdout'
petscopts['snes_rtol'] = '1e-4'


sys = System(nsi)
state = sys.x.duplicate()
time = numpy.linspace(0, .75, 100)

alpha = 0.5
for i in xrange(len(time)-1):
    print 'step: %d' % i
    sys.step(state, time[i], time[i+1], alpha)

dofmap = nsi.getDofMap()
stt = numpy.array(state)
sol = numpy.empty(nsi.getSize(), dtype=float)
dofmap.apply(stt, sol, time[-1])

sol.shape = (sol.size/2, 2)

newnod = xnod + sol;
out = file('xnod.out', 'w')
for xy in newnod:
    out.write('%f %f\n' % tuple(xy))
out.close()

out = file('icone.out', 'w')
for ijk in icone:
    out.write('%d %d %d\n' % tuple(ijk))
out.close()
