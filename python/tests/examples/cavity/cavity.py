if __name__ == '__main__':
    import sys, petsc4py
    petsc4py.init(sys.argv)
    del sys, petsc4py

import numpy
from petsc4py import PETSc
from pf4py.kernel import *
from cavmsh import Cavity

class Log:
    def __init__(self):
        self.total = 0.0
    def beg(self, msg=''):
        PETSc.Sys.Print('%s: start ...\n' % msg)
        self.tic = PETSc.Log.getTime()
    def end(self, msg=''):
        self.tac = PETSc.Log.getTime()
        elapsed = self.tac - self.tic
        self.total += elapsed
        PETSc.Sys.Print('%s: elapsed: %f, total: %f\n' % (msg, elapsed, self.total) )

log = Log()

opts = PETSc.Options()

D   = opts.getInt('D', 2)
N   = opts.getInt('N', 32)
rho = opts.getScalar('rho', 1.0)
Re  = opts.getScalar('Re',  1.0)
CFL = opts.getReal('CFL',   1.0)

SHAPE = (N, N, N)[:D]

xnod, icone, wall = Cavity(SHAPE)
nodedata = DTableS(xnod)
elemconn = DTableI(icone+1)
#del xnod, icone

nnod, ndim = nodedata.getShape()
ndof = ndim+1


domain = Domain(ndim, nnod, ndof)

domain.setNodedata(nodedata)

elemset = Elemset('nsi_tet_les_full', elemconn)
elemset.setOptions({
    'chunk_size'     : str(100000),
    'block_uploading': str(1),
    
    #'ndim'       : str(ndim),
    #'geometry'   : 'cartesian%dd' % ndim,
    #'npg'        : str(2**ndim),

    #'LES'    : str(1),
    #'C_smag' : str(0.18),

    'rho'        : str(rho),
    'viscosity'  : str(rho/Re),

    'use_lsic' : str(1),
    'temporal_stability_factor' : str(1),

    })

domain.addElemset(elemset)

# boundary conditions
# + pressure
#domain.setFixation(0, ndim, 0.0)            # p  = 0 (idof 3)

# + velocity
for k in wall:
    nods = wall[k]
    if k != 'top': 
        domain.setFixation(nods, 0, 0.0) # vx = 0 
    if ndim > 1:
        domain.setFixation(nods, 1, 0.0) # vy = 0 
    if ndim > 2:
        domain.setFixation(nods, 2, 0.0) # vz = 0
nods = wall['top']
#vlid =  - (xlid**elid + (xlid-1)**elid) + 1.0
vlid = 1
domain.setFixation(nods, 0, vlid) # vx = 1

domain.setUp()

steady = opts.getBool('steady', False)

if steady:
    from pf4py.solvers import SteadySolver as Solver
else:
    from pf4py.solvers import TransientSolver as Solver

if opts.getBool('schur', False):
    #solver = Solver(domain, mat_type='is')
    solver = Solver(domain)
    opts["ksp_type"] = "preonly"
    opts["pc_type"]   = "schur"
    opts["pc_schur_local_ccsize"] = 5000
    opts["pc_schur_print_stats"]  = 1

    opts["sub_ksp_type"] = "gmres"
    opts["sub_ksp_gmres_restart"] = 1000
    opts["sub_pc_type"]       = "bjacobi"
    #opts["sub_ksp_right_pc"]  = 1

    if 'ksp_monitor' in opts:
        opts['sub_ksp_monitor'] = opts['ksp_monitor']
        del opts['ksp_monitor']
        
elif opts.getBool('asm', False):
    solver = Solver(domain)
    ostr ="""
    -ksp_right_pc
    -pc_type asm
    #-pc_asm_overlap 10
    -ksp_gmres_restart 100
    """.split(); ostr = ' '.join(ostr)
    PETSc.Options.insertString(ostr)
else:
    solver = Solver(domain)
    ostr ="""
    -ksp_right_pc
    -ksp_gmres_restart 100
    """.split(); ostr = ' '.join(ostr)
    PETSc.Options.insertString(ostr)
    
    
if steady:
    snes = solver.snes
else:
    ts = solver.ts
    snes = ts.snes

snes.rtol  = 1e-6
snes.atol  = 1e-12
    
ksp  = snes.ksp

if 'ew' not in opts:
    ksp.rtol = 1e-3
else:
    #snes.setUseEW(True, version=1)
    #snes.setUseEW(True, version=2)
    snes.setUseEW(True, version=3)
    #snes.setParamsEW(rtol_0=0.5)
    #snes.setParamsEW(version=3,rtol_0=1e-3)
    #snes.setParamsEW(version=3,rtol_0=1e-2)
    

if steady:
    snes.setFromOptions()
    u0 = None
else:
    ts = solver.ts

    L  = 1.0
    h  = L/(N-1)
    U  = 1.0
    dT = CFL*h/abs(U)
    
    T0 = opts.getScalar('T0', 0.00)
    dT = opts.getScalar('dT', dT)
    nT = opts.getInt('nT', 5)
    T = T0 + nT*dT*1.001

    ts.setFromOptions()

    ts.setTime(T0)
    ts.setTimeStep(dT/CFL*2)
    ts.setDuration(T, 1)
    solver.run()

    u0=solver.solution

    ts.setTime(T0)
    ts.setTimeStep(dT)
    ts.setDuration(T, nT)


class Mon:
    def __init__(self):
        self.save = False
        self.count = 0
        self.rnorm  = []
        self.fnorm  = []
        self.index  = []

    def ts(self, ts, i, t, u):
        if i == 0: return
        self.save = True
        
    def snes(self, snes, i, fnorm):
        if not self.save : return
        self.fnorm.append(fnorm)
        self.index.append(self.count)
        
    def ksp(self, ksp, i, rnorm):
        if not self.save : return
        self.count += 1
        self.rnorm.append(rnorm)
        
mon = Mon()
if not steady:
    ts.setMonitor(mon.ts)
    snes.setMonitor(mon.snes)
    ksp.setMonitor(mon.ksp)
else:
    #mon.save = True
    #snes.setMonitor(mon.snes)
    #ksp.setMonitor(mon.ksp)
    pass

log.beg('solver.run()')
sol = solver.run(u0)
log.end('solver.run()')

lrh = mon.rnorm
lit = range(0, len(lrh))
nrh = mon.fnorm
nit = mon.index

PETSc.Sys.Print("nni:%d, nli:%d\n" %(len(nit),len(lit)))


sol = sol.copy()

x, y = numpy.mgrid[0:1:N*1j,0:1:N*1j]

vx = sol[:,0].reshape(N,N)
vy = sol[:,1].reshape(N,N)
u  = numpy.sqrt(vx**2 + vy**2)
p  = sol[:,2].reshape(N,N)

print p.max(), p.min(), p.sum()

vx1 = vx[:,N//2]

xx = yy = numpy.linspace(0,1,N)

#raise SystemExit

if PETSc.COMM_WORLD.getSize() > 1: raise SystemExit

try:
    from matplotlib import pylab
except ImportError:
    print "'matplotlib.pylab' not available"
    raise SystemExit

pylab.figure()
pylab.quiver(x,y,vx,vy,scale=5)
pylab.axis('equal')

pylab.figure()
pylab.contourf(x,y,u)
pylab.axis('equal')

pylab.figure()
pylab.contourf(x,y,p)
pylab.axis('equal')

pylab.show()
