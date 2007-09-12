import numpy

def mkCoord(shape, bbox):
    slices = tuple(slice(ax[0], ax[1], shape[i]*1j)
                   for i, ax in enumerate(bbox))
    grid = numpy.mgrid[slices]
    axes = [axis.ravel() for axis in grid]
    xnod = numpy.column_stack(axes)
    return xnod

def mkGrid(shape):
    nnod = numpy.product(shape)
    grid = numpy.arange(nnod)
    grid.shape = shape
    return grid

def mkQuad(grid):
    M, N = grid.shape
    econn =  numpy.zeros((M-1, N-1, 4), dtype=int)
    econn[..., 0] = grid[ 1:M,   0:N-1 ]
    econn[..., 1] = grid[ 1:M,   1:N   ]
    econn[..., 2] = grid[ 0:M-1, 1:N   ]
    econn[..., 3] = grid[ 0:M-1, 0:N-1 ]
    econn.shape = (econn.size/4, 4)
    return econn

def mkTri(grid):
    quad = mkQuad(grid)
    split = [1,2,3,
             3,4,1]
    split = [i-1 for i in split]
    tri = quad[:, split].copy()
    tri.shape = (tri.size/3, 3)
    return tri

def mkHexa(grid):
    N0, N1, N2 = grid.shape
    econn = numpy.zeros((N0-1, N1-1, N2-1, 8), dtype=int)
    econn[..., 0] = grid[ 0:N0-1 , 0:N1-1, 0:N2-1 ]
    econn[..., 1] = grid[ 0:N0-1 , 1:N1  , 0:N2-1 ]
    econn[..., 2] = grid[ 0:N0-1 , 1:N1  , 1:N2   ]
    econn[..., 3] = grid[ 0:N0-1 , 0:N1-1, 1:N2   ]
    econn[..., 4] = grid[ 1:N0   , 0:N1-1, 0:N2-1 ]
    econn[..., 5] = grid[ 1:N0   , 1:N1  , 0:N2-1 ]
    econn[..., 6] = grid[ 1:N0   , 1:N1  , 1:N2   ]
    econn[..., 7] = grid[ 1:N0   , 0:N1-1, 1:N2   ]
    econn.shape = (econn.size/8, 8)
    return econn

def mkTetra(grid):
    hexa = mkHexa(grid)
    split = [1,2,4,8,
             1,2,8,5,#1,2,5,8,
             2,3,4,8,
             2,3,8,7,#2,3,7,8,
             2,5,6,8,
             2,6,7,8,]
    split = [i-1 for i in split]
    tetra = hexa[:, split].copy()
    tetra.shape = (tetra.size/4, 4)
    return tetra


def SqCav(shape=(10, 10), bbox=((0, 1), (0, 1)),topo='quad'):
    grid = mkGrid(shape)
    xnod = mkCoord(shape, bbox)
    if topo == 'quad':
        econn = mkQuad(grid)
    elif topo == 'tri':
        econn = mkTri(grid)
    else:
        raise ValueError("bad 'topo' %s" % topo)
    wall = {'left'   : grid[ 0,  :],
            'right'  : grid[-1,  :],
            'bottom' : grid[ :,  0],
            'top'    : grid[ :, -1],}
    return xnod, econn, wall

def CubCav(shape=(5, 5, 5), bbox=((0, 1), (0, 1), (0, 1)),topo='hexa'):
    grid = mkGrid(shape)
    xnod = mkCoord(shape, bbox)
    econn = mkHexa(grid)
    if topo == 'hexa':
        econn = mkHexa(grid)
    elif topo == 'tetra':
        econn = mkTetra(grid)
    else:
        raise ValueError("bad 'topo' %s" % topo)
    wall = {'left'   : grid[  0,    ...],
            'right'  : grid[ -1,    ...],
            'front'  : grid[  :,  0,  :],
            'back'   : grid[  :, -1,  :],
            'bottom' : grid[...,      0],
            'top'    : grid[...,     -1],}
    return xnod, econn, wall

def Cavity(shape, *targs, **kargs):
    if len(shape) == 2:
        return SqCav(shape, *targs, **kargs)
    elif len(shape) == 3:
        return CubCav(shape, *targs, **kargs)
    else:
        raise ValueError('invalid shape')
