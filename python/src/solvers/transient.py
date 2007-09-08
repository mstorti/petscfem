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


class TEqSys:

    ALPHA     = 1.0
    STRUCTURE = PETSc.Mat.Structure.SAME_NONZERO_PATTERN

    def __init__(self, alpha=None):
        if alpha is None:
            alpha = TEqSys.ALPHA
        alpha = float(alpha)
        assert alpha >= 0.0
        assert alpha <= 1.0
        self.alpha = alpha
        self._skip_res = False
    
    def computeResidual(self, ts, t, u, R):
        t0 = t - ts.getTimeStep()
        u0 = ts.getSolution()
        domain = ts.getAppCtx()
        R.zeroEntries()
        domain.assemble(t, u, t0, u0, R, None, self.alpha)
        
    def computeJacobian(self, ts, t, u, J, P):
        t0 = t - ts.getTimeStep()
        u0 = ts.getSolution()
        domain = ts.getAppCtx()
        P.zeroEntries()
        domain.assemble(t, u, t0, u0, None, P, self.alpha)
        if J != P: J.assemble()
        return self.STRUCTURE


class TransientSolver:

    """
    Transient Solver, based in PETSc TS timestepper solver.
    """
    
    def __init__(self, domain, alpha=None, mat_type=''):
        # init  members
        self.eqsys = TEqSys(alpha)
        self.solution = PETSc.Vec()
        self.state = PETSc.Vec()
        self.residual = PETSc.Vec()
        self.jacobian = PETSc.Mat()
        self.ts = PETSc.TS()
        # create sub objects
        domain.allocateSolution(self.solution)
        domain.allocateState(self.state)
        domain.allocateResidual(self.residual)
        domain.allocateJacobian(self.jacobian, mat_type)
        self.ts.create(domain.getComm())
        self.ts.setProblemType('NONLINEAR')
        self.ts.setType('user')
        self.ts.setSolution(self.state)
        self.ts.setFunction(self.eqsys.computeResidual, self.residual)
        self.ts.setJacobian(self.eqsys.computeJacobian, self.jacobian)
        self.ts.setAppCtx(domain)
        
    def run(self, u0=None):
        from numpy import asarray
        domain = self.ts.getAppCtx()
        # build initial solution
        if u0 is not None:
            u0 = asarray(u0, dtype=PETSc.Scalar)
            self.solution.setArray(u0)
            time = self.ts.getTime()
            domain.buildState(time, self.solution, self.state)
        else:
            self.state.set(0)
        # solve and build solution
        self.ts.solve(self.state)
        time = self.ts.getTime()
        domain.buildSolution(time, self.state, self.solution)
        # return solution as an array
        from numpy import asarray
        sol = asarray(self.solution)
        nnod = domain.getNNod()
        ndof = domain.getNDof()
        sol.shape = (nnod, ndof)
        return sol
