from petscfem4py import Nodeset

opts = dict(a='0', b='1', c='2')

n1 = Nodeset()
assert n1.getRefCount() == 1
assert n1.getDataSize() == (0, 0)
assert n1.getData().shape == (0, 0)
assert n1.getOptions() == { }
n1.setOptions(opts)
n1.setData([[0,1,2] for i in xrange(10)])
assert n1.getOptions() == opts
assert n1.getDataSize() == (10, 3)
assert n1.getData().shape == (10, 3)


n2 = Nodeset(3, [[0,1,2] for i in xrange(10)])
assert n2.getRefCount() == 1
assert n2.getOptions() == { }
n2.setOptions(opts)
assert n2.getOptions() == opts
assert n2.getDataSize() == (10, 3)
assert n2.getData().shape == (10, 3)


n3 = Nodeset(n2)
assert n3 != n2
assert n3.getRefCount() == 1
assert n3.getOptions()  == n2.getOptions()
assert n3.getDataSize() == n2.getDataSize()
assert (n3.getData()    == n2.getData()).all()


n4 = Nodeset()
n4.setData(n3.getData())
n4.setOptions(n3.getOptions())
assert n4.getRefCount() == 1
assert n4.getOptions()  == n3.getOptions()
assert n4.getDataSize() == n3.getDataSize()
assert (n4.getData()    == n3.getData()).all()


## n1.clear()
## n2.clear()
## n3.clear()
## n4.clear()
