from petscfem4py.petscfem import *

dofset = Dofset(6, 3)

fixation = ([ 1,   3,   5],
            [ 0,   1,   2],
            [77., 88., 99.])

constraint = ([  0,    2,     4],
              [  0,    1,     2],
              [250., 500., 1000.])

dofset.addFixations(*fixation)
dofset.addConstraints(*constraint)
