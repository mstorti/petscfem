# Author:  Lisandro Dalcin
# Contact: dalcinl@users.sourceforge.net
# License: Public Domain
# $Id$

"""
PETSc-FEM for Python
"""

__author__    = 'Lisandro Dalcin <dalcinl@intec.unl.edu.ar>'
__credits__   = "Mario Storti <mstorti@intec.unl.edu.ar>"

__date__      = '$Date$'
__version__   = '$Version$'
__revision__  = '$Revision$'

__docformat__ = 'reST'


from petsc4py import PETSc


class SEqSys:

    STRUCTURE = PETSc.Mat.Structure.SAME_NONZERO_PATTERN

    def __init__(self):
        self._skip_res = False

    def computeResidual(self, snes, x, R):
        R.zeroEntries()
        if self._skip_res:
            self._skip_res = False
            return
        J, P, _ = snes.getJacobian()
        mffd = J.getType() == PETSc.Mat.Type.MFFD
        domain = snes.getApplicationContext()
        if mffd:
            domain.assemble(0.0, x, R, None)
        else:
            J.zeroEntries()
            domain.assemble(0.0, x, R, J)
        
    def computeJacobian(self, snes, x, J, P):
        if snes.its == snes.max_it - 1:
            self._skip_res = True
        J, P, _ = snes.getJacobian()
        mffd = J.getType() == PETSc.Mat.Type.MFFD
        domain = snes.getApplicationContext()
        if mffd:
            P.zeroEntries()
            domain.assemble(0.0, x, None, P)
            J.assemble()
        return self.STRUCTURE


class SteadySolver:

    """
    Steady-State Solver, based in PETSc SNES nonlinear solver.
    """
    
    def __init__(self, domain, mat_type=''):
        # init  members
        self.eqsys  = SEqSys()
        self.solution = PETSc.Vec()
        self.state = PETSc.Vec()
        self.residual = PETSc.Vec()
        self.jacobian = PETSc.Mat()
        self.snes = PETSc.SNES()
        # create sub objects
        domain.allocateSolution(self.solution)
        domain.allocateState(self.state)
        domain.allocateResidual(self.residual)
        domain.allocateJacobian(self.jacobian,mat_type)
        self.snes.create(domain.getComm())
        #self.snes.setType('ls')
        #snesls = PETSc.SNESLS(self.snes)
        #snesls.setLineSearch('basic')
        self.snes.setFunction(self.eqsys.computeResidual, self.residual)
        self.snes.setJacobian(self.eqsys.computeJacobian, self.jacobian)
        self.snes.setApplicationContext(domain)
        
    def run(self, u0=None):
        from numpy import asarray
        domain = self.snes.getApplicationContext()
        # build initial solution
        if u0 is not None:
            u0 = asarray(u0, dtype=PETSc.Scalar)
            self.solution.setArray(u0)
            time = 0.0
            domain.buildState(time, self.solution, self.state)
        else:
            self.state.set(0)
        # solve nonlinear system
        self.snes.solve(self.state)
        # build final solution
        domain.buildSolution(0.0, self.state, self.solution)
        # return the solution as a numpy array
        # with appropriate shape
        asol = asarray(self.solution)
        nnod = domain.getNNod()
        ndof = domain.getNDof()
        asol.shape = (nnod, ndof)
        return asol
