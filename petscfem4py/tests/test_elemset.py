from petscfem4py.petscfem import *

icone = [[0, 1, 2, 3]]

elemset = Elemset('nsi_tet_les_fm2')

elemset.setConnectivity(icone)

icone = elemset.getConnectivity()

print icone
