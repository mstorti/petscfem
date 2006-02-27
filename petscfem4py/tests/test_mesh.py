from petscfem4py.petscfem import *

xnod  = [[0,0],
         [1,0],
         [1,1],
         [0,1]]

icone = [[0, 1, 2, 3]]

nodedata = Nodedata(xnod)

elemset1 = Elemset('nsi_tet_les_fm2', 'e0')
elemset1.setConnectivity(icone)

elemset2 = Elemset('nsi_tet_les_comp', 'e1')
elemset2.setConnectivity(icone)

elemset3 = Elemset('nsi_tet_les_ther', 'e2')
elemset3.setConnectivity(icone)

mesh = Mesh()
mesh.setNodedata(nodedata)
for e in [elemset1, elemset2, elemset3]:
    mesh.addElemset(e)

e0 = mesh.getElemset(0)
e1 = mesh.getElemset(1)
e2 = mesh.getElemset(2)

assert mesh.hasElemset('e0')
assert mesh.hasElemset('e1')
assert mesh.hasElemset('e2')

e0 = mesh.findElemset('e0')
e1 = mesh.findElemset('e1')
e2 = mesh.findElemset('e2')
