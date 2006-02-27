from petscfem4py.petscfem import *

dm = DofMap(6, 3)

fixa = ([ 1,   3,   5],
        [ 0,   1,   2],
        [77., 88., 99.])

dm.addFixations(fixa)


constraint = ([  0,    2,     4],
              [  0,    1,     2],
              [250., 500., 1000.])

dm.addConstraints(constraint)

# dm.view()
