from petscfem4py import *

icone = [[0, 1, 2, 3]]

elemset = Elemset('nsi_tet_les_fm2', icone)
#elemset.setName('navier-stokes')
#elemset.setNDof(3)

## for i in range(4):
##     elemset.setOption('key%d'%i, 'value%d'%i)
## for i in range(4):
##     val = elemset.getOption('key%d'%i)
##     assert val == 'value%d'%i, 'invalid value' 

## assert elemset.getType() == 'nsi_tet_les_fm2'
## assert elemset.getName() == 'navier-stokes'

## print elemset.getData()
## elemset.view()
