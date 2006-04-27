import numpy
import petsc.PETSc as PETSc
from petscfem4py import *

ndim = 2
ndof = 3

xnod = """\
0.00000 0.00000
0.00000 0.10000
0.10000 0.00000
0.10000 0.10000
0.20000 0.00000
0.20000 0.10000
0.30000 0.00000
0.30000 0.10000
0.40000 0.00000
0.40000 0.10000
0.50000 0.00000
0.50000 0.10000
0.60000 0.00000
0.60000 0.10000
0.70000 0.00000
0.70000 0.10000
0.80000 0.00000
0.80000 0.10000
0.90000 0.00000
0.90000 0.10000
1.00000 0.00000
1.00000 0.10000
""".split()
xnod = numpy.array([float(i) for i in xnod])
xnod.shape = (xnod.size/ndim, ndim)
nodeset = Nodeset(ndim, xnod)

icone = """\
 0  2  3  1
 2  4  5  3
 4  6  7  5
 6  8  9  7
 8 10 11  9
10 12 13 11
12 14 15 13
14 16 17 15
16 18 19 17
18 20 21 19
""".split()
icone = numpy.array([int(i) for i in icone])
icone.shape = (icone.size/4, 4)
elemset = Elemset('nsi_tet_les_fm2', icone)
elemset.setName('fluid')
opts ={'geometry'  : 'cartesian2d',
       'ndim'      : str(ndim),
       'npg'       : str(4),
       'viscosity' : str(100),
       'weak_form' : str(0),
       }
elemset.setOptions(opts)


dofset = Dofset(nodeset.getSize(), ndof)

# velocity and pressure set to zero
for node in xrange(20, 22):
    for field in xrange(0, ndof):
        dofset.addFixations(node, field, 0.0)

# time dependet boundary conditions
class Sin(Amplitude):
    def __init__(self, A, w):
        super(Sin, self).__init__()
        self.A = A
        self.w = w
    def __call__(self, node, field, time):
        #print node, field, time
        #print time, self.A * numpy.sin(self.w*time)
        return self.A * numpy.sin(self.w*time)

amp = Sin(1.0, 2*numpy.pi)

for node in xrange(0, 2):
    dofset.addFixations(node, 0, 0)
    dofset.addFixations(node, 1, 1, amp)

# periodic boundary conditions from y=0 to y=0.1.
nodes  = [0,  1]
fields = [2,  2]
coeffs = [1, -1]
dofset.addConstraints(nodes, fields, coeffs)
for node in xrange(2, 20, 2):
    nodes  = [node, node+1]
    for field in xrange(0, ndof):
        fields = [field, field]
        dofset.addConstraints(nodes, fields, coeffs)

   
nsi = NavierStokes(nodeset, [elemset], dofset)

class System:

    def __init__(self, problem):
        self.problem = problem

        COMM  = problem.getComm()
        comm  = COMM.comm
        ndofs = problem.getDofSizes()
        if COMM.size == 1:
            J  = PETSc.MatSeqAIJ(ndofs, nz=9*4, comm=comm)
        else:
            J  = PETSc.MatMPIAIJ(ndofs, d_nz=9*4, comm=comm)
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
if nsi.comm.size==1:
    petscopts['ksp_type']     = 'preonly'
    petscopts['pc_type']      = 'lu'
petscopts['snes_monitor'] = 'stdout'
petscopts['ksp_monitor']  = 'stdout'
#petscopts['snes_view']    = 'stdout'


sys = System(nsi)
state = sys.x.duplicate()
time = numpy.linspace(0, 2.5, 11)

alpha = 0.5
for i in xrange(len(time)-1):
    sys.step(state, time[i], time[i+1], alpha)
    print state.norm()

dofmap = nsi.getDofMap()
stt = numpy.array(state)
sol = numpy.empty(nsi.getSize(), dtype=float)
dofmap.apply(stt, sol, time[-1])

out = numpy.fromfile('oscplate.out', dtype=float, sep=' ')

sol.shape = (sol.size/ndof, ndof)
out.shape = (out.size/ndof, ndof)
