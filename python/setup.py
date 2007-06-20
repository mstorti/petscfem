#!/bin/env python

# $id$

name    = 'petscfem4py'

version = '0.1.0'

classifiers = """
License :: Public Domain
Operating System :: POSIX
Intended Audience :: Developers
Intended Audience :: Science/Research
Programming Language :: C
Programming Language :: C++
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries :: Python Modules
"""

keywords = """
scientific computing
parallel computing
finite elements method
"""

metadata = {
    'name'             : name,
    'version'          : version,
    'author'           : 'Lisandro Dalcin',
    'author_email'     : 'dalcinl@users.sourceforge.net',
    'url'              : 'http://www.cimec.org.ar/python',
    'classifiers'      : [c for c in classifiers.split('\n') if c],
    'keywords'         : [k for k in keywords.split('\n')    if k],
    'license'          : 'Public Domain',
    'platforms'        : ['POSIX'],
    'maintainer'       : 'Lisandro Dalcin',
    'maintainer_email' : 'dalcinl@users.sourceforge.net',
    }

# --------------------------------------------------------------------
# Basic Configuration
# --------------------------------------------------------------------
from os import environ as env
from posixpath import expanduser, abspath, join

MPI_DIR      = env.get('MPI_DIR')    or '/usr/local/mpich2/1.0.5'
PETSC_DIR    = env.get('PETSC_DIR')  or '/usr/local/petsc/2.3.3'
PETSC_ARCH   = env.get('PETSC_ARCH') or 'linux-gnu'
PETSCFEM_DIR = env.get('PETSCFEM_DIR') or abspath('..')
SOFT         = abspath(join(PETSCFEM_DIR, '..', 'SOFT'))

GLIB1    = '/usr/lib/glib/include'
GLIB2    = '/usr/include/glib-1.2'
LIBRETTO = join(SOFT, 'libretto-2.1')
NEWMAT   = join(SOFT, 'newmat-1.0')
MESCHACH = join(SOFT, 'meschach-1.2')
METIS    = join(SOFT, 'metis-4.0')
ANN      = join(SOFT, 'ann-1.1.1')


BOPT = env.get('BOPT') or 'g'
if PETSC_ARCH.endswith('O_c++') or  PETSC_ARCH.endswith('O'):
    BOPT = 'O'
if PETSC_ARCH.endswith('g_c++') or PETSC_ARCH.endswith('g'):
    BOPT = 'g'

#MKL_DIR = ['/opt/intel/mkl/9.0/lib/32']
#MKL_LIB = ['mkl_lapack', 'mkl_def', 'guide', 'vml', 'pthread']
MKL_DIR = []
MKL_LIB = []

# --------------------------------------------------------------------
# Extension modules
# --------------------------------------------------------------------

config = {
    'define_macros': [('MPICH_SKIP_MPICXX', None),
                      ('USE_ANN', None),
                       ],
    'include_dirs' : [join(MPI_DIR, 'include'),
                      join(PETSC_DIR, 'bmake', PETSC_ARCH),
                      join(PETSC_DIR, 'include'),
                      join(LIBRETTO, 'include'),
                      join(NEWMAT, 'src'),
                      join(ANN, 'include'),
                      GLIB1, GLIB2,
                      MESCHACH, METIS,
                      ],
    
    'libraries'    : ['mpich',
                      'petscksp', 'petscdm', 'petscmat', 'petscvec', 'petsc',
                      'glib',
                      'newmat',
                      'ibretto',
                      'meschach',
                      'metis',
                      'ANN',
                      ],
    
    'library_dirs' : [join(MPI_DIR, 'lib'),
                      join(PETSC_DIR, 'lib'),
                      join(PETSC_DIR, 'lib', PETSC_ARCH),
                      join(LIBRETTO, 'lib'),
                      join(ANN, 'lib'),
                      join(NEWMAT, 'src'),
                      MESCHACH, METIS,
                      ],
    'runtime_library_dirs' : []
    }

def ext_modules(Extension):
    from os.path import abspath, join

    PETSCFEM_MACROS = config['define_macros']
    PETSCFEM_DIR = abspath('..')
    PETSCFEM_INCLUDE = [PETSCFEM_DIR, join(PETSCFEM_DIR, 'src')] \
                       + config['include_dirs']

    if BOPT in ('g', 'g_c++'):
        PETSCFEM_LIBRARY = ['ns_g', 'petscfem_g'] * 2
    elif BOPT in ('O',  'O_c++'):
        PETSCFEM_LIBRARY = ['ns_O', 'petscfem_O'] * 2
    else:
        raise SystemExit('invalid BOPT: %s' % BOPT)
    PETSCFEM_LIBRARY +=  config['libraries']
    PETSCFEM_LIBDIR = [join(PETSCFEM_DIR, 'src'),
                       join(PETSCFEM_DIR, 'applications/ns') ] \
                            + config['library_dirs']
    PETSCFEM_RT_LIBDIR = PETSCFEM_LIBDIR \
                          + config['runtime_library_dirs']
    
    from glob import glob
    sources = glob('src/kernel/core.i')
    sources += glob('src/kernel/*.cpp')
    try: 
        sources.remove('src/kernel/core_wrap.cpp')
    except ValueError:
        pass
    depends = glob('src/kernel/*.i') + glob('src/kernel/*.h')

    PETSCFEM_LIBRARY   += MKL_LIB
    PETSCFEM_LIBDIR    += MKL_DIR
    PETSCFEM_RT_LIBDIR += MKL_DIR

    petscfem = Extension('pf4py.kernel._core',
                         sources=sources,
                         depends=depends,
                         define_macros        = PETSCFEM_MACROS,
                         include_dirs         = PETSCFEM_INCLUDE,
                         libraries            = PETSCFEM_LIBRARY,
                         library_dirs         = PETSCFEM_LIBDIR,
                         runtime_library_dirs = PETSCFEM_RT_LIBDIR,
                         language='c++',
                       )
    return [petscfem]

def c_libraries():
    from os.path import abspath, join

    PETSCFEM_DIR = abspath('..')
    PETSCFEM_MACROS = config['define_macros']
    PETSCFEM_INCLUDE_DIR = [PETSCFEM_DIR, join(PETSCFEM_DIR, 'src')] \
                           + config['include_dirs']
    from glob import glob
    sources = glob('src/kernel/*.cpp')
    try: 
        sources.remove('src/kernel/core_wrap.cpp')
    except ValueError:
        pass
    depends = glob('src/kernel/*.h')

    return [('kernel', dict(sources      = sources,
                            depends      = depends,
                            macros       = PETSCFEM_MACROS,
                            include_dirs = PETSCFEM_INCLUDE_DIR))]

# --------------------------------------------------------------------
# Setup
# --------------------------------------------------------------------

def setup():

    from numpy.distutils.core import setup
    from numpy.distutils.core import Extension
    from numpy.distutils.command.build_src import build_src as _build_src

    class build_src(_build_src):
        def finalize_options(self):
            _build_src.finalize_options(self)
            self.inplace = True
            try:
                swig_opts = self.swig_opts
            except AttributeError:
                swig_opts = self.swigflags
            for flag in ('-O',
                         #'-nodefaultctor',
                         #'-nodefaultdtor',
                         #'-keyword',
                         '-nofastproxy',
                         '-noh',
                         '-nodirvtable', # bug in swig 1.3.29
                         ):
                swig_opts.append(flag)
            
            
    setup(packages     = ['pf4py',
                          'pf4py.kernel',
                          'pf4py.solvers',
                          'pf4py.tools',
                          ],
    	  package_dir  = {'pf4py' : 'src'},
          ext_modules  = ext_modules(Extension),
          #libraries   = c_libraries(),
          cmdclass     = {'build_src' : build_src},
          **metadata)

    
# --------------------------------------------------------------------

if __name__ == '__main__':
    
    from distutils import sysconfig
    cvars = sysconfig.get_config_vars()
    cflags = cvars['OPT'].split()
    for flag in ('-Wall',
                 '-Wstrict-prototypes',):
        try:
            cflags.remove(flag)
        except ValueError:
            pass
    if BOPT in ('g',  'g_c++'):
        try:
            cflags.remove('-DNDEBUG')
        except ValueError:
            pass
    elif BOPT in ('O',  'O_c++'):
        try:
            cflags.remove('-g')
        except ValueError:
            pass
        cflags.append('-funroll-loops')
    else:
        raise SystemExit('invalid BOPT: %s' % BOPT)
    
    cvars['OPT'] = str.join(' ', cflags) 

    setup()

# --------------------------------------------------------------------
