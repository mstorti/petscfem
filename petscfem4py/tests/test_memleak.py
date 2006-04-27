import numpy
from petscfem4py import *


def memusage():
    import os
    _proc_pid_stat = '/proc/%s/stat' % os.getpid()
    try:
        f=open(_proc_pid_stat, 'r')
        l = f.readline().split(' ')
        f.close()
        return int(l[22])
    except Exception:
        return 0
        

print memusage()
print
print

print memusage()
obj = Nodeset()
print memusage()
opt = dict((str(i), str(i)) for i in xrange(10000))
print memusage()
for i in xrange(5):
    obj.setOptions(opt)
    obj.clear()
    print memusage()
del obj, opt
print memusage()
print
print


print memusage()
xnod = numpy.ones((1000000,3), dtype=float)
print memusage()
nodeset = Nodeset(3, xnod)
print memusage()
for i in xrange(5):
    nodeset.setData(xnod)
    nodeset.clear()
    print memusage()
for i in xrange(5):
    nn = Nodeset(nodeset)
    del nn
    print memusage()
del xnod, nodeset
print memusage()
print
print
