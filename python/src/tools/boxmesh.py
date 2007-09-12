import numpy

def mkGrid(shape):
    nnod = numpy.product(shape)
    grid = numpy.arange(nnod)
    grid.shape = shape
    return grid

def mkCoord(shape, bbox):
    slices = tuple(slice(ax[0], ax[1], shape[i]*1j)
                   for i, ax in enumerate(bbox))
    grid = numpy.mgrid[slices]
    axes = [axis.ravel() for axis in grid]
    xnod = numpy.column_stack(axes)
    return xnod

def mkQuad(grid):
    M, N = grid.shape
    quads =  numpy.zeros((M-1, N-1, 4), dtype=int)
    quads[..., 0] = grid[ 1:M,   0:N-1 ]
    quads[..., 1] = grid[ 1:M,   1:N   ]
    quads[..., 2] = grid[ 0:M-1, 1:N   ]
    quads[..., 3] = grid[ 0:M-1, 0:N-1 ]
    quads.shape = (-1, 4)
    return quads

def mkTri(grid):
    quads = mkQuad(grid)
    split = [1,2,3,
             3,4,1]
    split = [i-1 for i in split]
    tri = quads[:, split].copy()
    tri.shape = (-1, 3)
    return tri

def mkHexa(grid):
    N0, N1, N2 = grid.shape
    hexas = numpy.zeros((N0-1, N1-1, N2-1, 8), dtype=int)
    hexas[..., 0] = grid[ 0:N0-1 , 0:N1-1, 0:N2-1 ]
    hexas[..., 1] = grid[ 0:N0-1 , 1:N1  , 0:N2-1 ]
    hexas[..., 2] = grid[ 0:N0-1 , 1:N1  , 1:N2   ]
    hexas[..., 3] = grid[ 0:N0-1 , 0:N1-1, 1:N2   ]
    hexas[..., 4] = grid[ 1:N0   , 0:N1-1, 0:N2-1 ]
    hexas[..., 5] = grid[ 1:N0   , 1:N1  , 0:N2-1 ]
    hexas[..., 6] = grid[ 1:N0   , 1:N1  , 1:N2   ]
    hexas[..., 7] = grid[ 1:N0   , 0:N1-1, 1:N2   ]
    hexas.shape = (-1, 8)
    return hexas

def mkTetra(grid):
    hexas = mkHexa(grid)
    split = [1,2,4,8,
             1,2,8,5,
             2,3,4,8,
             2,3,8,7,
             2,5,6,8,
             2,6,7,8,]
    split = [i-1 for i in split]
    tetras = hexas[:, split].copy()
    tetras.shape = (-1, 4)
    return tetras

def BoxMesh2D(shape, bbox=((0, 1), (0, 1)), topology='quad'):
    assert len(shape) == 2
    assert len(bbox)  == 2
    assert topology in ('quad', 'tri')
    grid = mkGrid(shape)
    xnod = mkCoord(shape, bbox)
    if topology == 'quad':
        icone = mkQuad(grid)
    elif topology == 'tri':
        icone = mkTri(grid)
    wall = {'left'   : grid[ 0,  :],
            'right'  : grid[-1,  :],
            'bottom' : grid[ :,  0],
            'top'    : grid[ :, -1],}
    return xnod, icone, wall

def BoxMesh3D(shape, bbox=((0, 1), (0, 1), (0, 1)), topology='hexa'):
    assert len(shape) == 3
    assert len(bbox)  == 3
    assert topology in ('hexa', 'tetra')
    grid = mkGrid(shape)
    xnod = mkCoord(shape, bbox)
    if topology == 'hexa':
        icone = mkHexa(grid)
    elif topology == 'tetra':
        icone = mkTetra(grid)
    wall = {'left'   : grid[  0,    ...],
            'right'  : grid[ -1,    ...],
            'front'  : grid[  :,  0,  :],
            'back'   : grid[  :, -1,  :],
            'bottom' : grid[...,      0],
            'top'    : grid[...,     -1],}
    return xnod, icone, wall

def BoxMesh(shape, **kargs):
    shape = tuple(shape)
    ndim = len(shape)
    if ndim not in (2, 3):
        raise ValueError('invalid dimension %d, only 2 or 3' % ndim)
    if ndim == 2:
        return BoxMesh2D(shape, **kargs)
    elif ndim == 3:
        return BoxMesh3D(shape, **kargs)



def _test_2d(N):
    x, e, w = BoxMesh([N,N])
    BoxMesh([N,N], bbox=[[-1,1], [-1,1]])
    BoxMesh([N,N], topology='quad')
    BoxMesh([N,N], topology='tri')
    return BoxMesh([N,N])

def _test_3d(N):
    BoxMesh([N,N,N])
    BoxMesh([N,N,N], bbox=[[-1,1], [-1,1], [-1,1]])
    BoxMesh([N,N,N], topology='hexa')
    BoxMesh([N,N,N], topology='tetra')
    return BoxMesh([N,N,N])

def _test(D, N):
    if D==2: return _test_2d(N)
    if D==3: return _test_3d(N)

if __name__ == '__main__':
    for D in (2, 3):
        for N in (1, 2, 3, 4):
            x, e, w = _test(D, N);
            assert x.shape == (   N**D,    D   )
            assert e.shape == ( (N-1)**D, 2**D )
