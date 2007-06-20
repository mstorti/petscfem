if __name__ == '__main__':
    import sys, petsc4py
    petsc4py.init(sys.argv)
    del sys, petsc4py

import numpy
from petsc4py import PETSc
from pf4py.kernel import *


xnod  = numpy.zeros([50, 3])

nnod, ndim = xnod.shape
ndof = ndim+1

nodedata = DTableS(xnod)
domain = Domain(ndim, nnod, ndof)
domain.setNodedata(nodedata)

field = DTableS(xnod)
mesh = domain.getMesh()
mesh.setField('field1', field)
mesh.setField('field2', field)
mesh.setField('field3', field)
mesh.setField('field4', field)

field1 = mesh.getField('field1')
field2 = mesh.getField('field2')
field3 = mesh.getField('field3')
field4 = mesh.getField('field4')

print field.__refcount__()

del field1, field2, field3, field4

print field.__refcount__()
