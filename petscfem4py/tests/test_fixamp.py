import numpy
from petsc import PETSc
from petscfem4py import *

nodeset = Nodeset(2, [[0,1] for i in xrange(4)])
elemset = Elemset('nsi_tet_les_fm2', [[0,1,2,3]])
options = """\
geometry   cartesian2d
ndim       2
npg        4
viscosity  0.001
weak_form  1"""
options = dict(line.split()
               for line in options.split('\n')
               if not line.startswith('#'))
elemset.setOptions(options)


mesh = Mesh(nodeset, [elemset])

dofset = Dofset(nodeset.getSize(), 3)

class MyAmp(Amplitude):
    def __call__(self, node, field, time):
        return node + field

amp = MyAmp()

dofset.addFixations(0, 0, 1)
dofset.addFixations(1, 0, 1)
dofset.addFixations(2, 0, 1, amp)
dofset.addFixations(3, 0, 1, amp)

problem = Problem(mesh, dofset)
dofmap = problem.getDofMap()

nnod = nodeset.getSize()
ndim = nodeset.getDim()
ndof = dofmap.getNDof();

sizes = problem.getDofSizes()
stt = PETSc.VecSeq(sizes)
sol = PETSc.VecSeq(nnod*ndof)

problem.buildSolution(stt, sol)
