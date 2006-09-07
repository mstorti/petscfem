from petscfem4py import Nodeset

opts = dict(a='0', b='1', c='2')

n1 = Nodeset(0, 3)
assert n1.__refcount__() == 1
assert n1.getSizes() == (0, 3, 0)
assert n1.getArray().shape == (0, 3)
assert n1.options == { }
n1.options = opts
n1.setArray([[0,1,2] for i in xrange(10)])
assert n1.options == opts
assert n1.getSizes() == (10, 3, 0)
assert n1.getArray().shape == (10, 3)


n2 = Nodeset([[0,1,2] for i in xrange(10)])
assert n2.__refcount__() == 1
assert n2.options == { }
n2.options = opts
assert n2.options == opts
assert n2.getSizes() == (10, 3, 0)
assert n2.getArray().shape == (10, 3)


n3 = Nodeset(n2)
assert n3 != n2
assert n3.__refcount__() == 1
assert n3.options  == n2.options
assert n3.getSizes() == n2.getSizes()
assert (n3.getArray() == n2.getArray()).all()


n4 = Nodeset(*n3.getSizes())
n4.setArray(n3.getArray())
n4.options = n3.options
assert n4.__refcount__() == 1
assert n4.options == n3.options
assert n4.getSizes() == n3.getSizes()
assert (n4.getArray() == n3.getArray()).all()


## n1.clear()
## n2.clear()
## n3.clear()
## n4.clear()
