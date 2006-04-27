from petscfem4py.petscfem import *

xnod  = [[0,0], [1,0], [1,1], [0,1]]
icone = [[0, 1, 2, 3] for i in xrange(4)]

nodedata = Nodeset(2, xnod)
elemset0 = Elemset('nsi_tet_les_fm2',  icone)
elemset1 = Elemset('nsi_tet_les_comp', icone)
elemset2 = Elemset('nsi_tet_les_ther', icone)

# default constructor
mesh1 = Mesh()
mesh1.setNodeset(nodedata)
for e in [elemset0, elemset1, elemset2]:
    mesh1.addElemset(e)
assert nodedata == mesh1.getNodeset()
assert elemset0 == mesh1.getElemset(0)
assert elemset1 == mesh1.getElemset(1)
assert elemset2 == mesh1.getElemset(2)

# copy constructor
mesh2 = Mesh(mesh1)
assert mesh2 != mesh1
assert nodedata == mesh1.getNodeset()
assert elemset0 == mesh1.getElemset(0)
assert elemset1 == mesh1.getElemset(1)
assert elemset2 == mesh1.getElemset(2)

# constructor from nodedata and elemset list
mesh3 = Mesh(nodedata, [elemset0, elemset1, elemset2])
assert nodedata == mesh1.getNodeset()
assert elemset0 == mesh1.getElemset(0)
assert elemset1 == mesh1.getElemset(1)
assert elemset2 == mesh1.getElemset(2)

mesh3.setNodeset(nodedata)
mesh3.setElemset(0, elemset0)
mesh3.setElemset(1, elemset1)
mesh3.setElemset(2, elemset2)

assert nodedata.getRefCount() == 4
assert elemset0.getRefCount() == 4
assert elemset1.getRefCount() == 4
assert elemset2.getRefCount() == 4

mesh1.clear()
mesh2.clear()
mesh3.clear()

assert nodedata.getRefCount() == 1
assert elemset0.getRefCount() == 1
assert elemset1.getRefCount() == 1
assert elemset2.getRefCount() == 1
