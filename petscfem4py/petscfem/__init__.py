# Author:  Lisandro Dalcin
# Contact: dalcinl@users.sourceforge.net
# License: GPL
# $Id: __init__.py,v 1.1.2.5 2006/05/24 21:12:21 dalcinl Exp $

"""
PETSc-FEM for Python
"""

__author__    = 'Lisandro Dalcin <dalcinl@intec.unl.edu.ar>'
__credits__   = "Mario Storti <mstorti@intec.unl.edu.ar>"

__date__      = '$Date: 2006/05/24 21:12:21 $'
__version__   = '$Version$'
__revision__  = '$Revision: 1.1.2.5 $'

__docformat__ = 'reST'


__all__ = ['OPTIONS',
           'Comm',
           'Object',
	   'Nodeset',
	   'Elemset',
	   'Mesh',
	   'Amplitude',
           'Dofset',
	   'DofMap',
           'Domain',

	   'Problem',
	   'NavierStokes',

           'Application',
           'NvrStks',
           ]

from petscfem import OPTIONS
from petscfem import Comm
from petscfem import Object
from petscfem import Nodeset
from petscfem import Elemset
from petscfem import Mesh
from petscfem import DofMap
from petscfem import Amplitude
from petscfem import Dofset
from petscfem import Domain

from petscfem import Problem
from petscfem import NavierStokes

from petscfem import Application
from petscfem import NvrStks
