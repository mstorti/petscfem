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
    'author_email'     : 'dalcinl@gmail.com',
    'url'              : 'http://www.cimec.org.ar/petscfem',
    'classifiers'      : [c for c in classifiers.split('\n') if c],
    'keywords'         : [k for k in keywords.split('\n')    if k],
    'license'          : 'Public Domain',
    'platforms'        : ['POSIX'],
    'maintainer'       : 'Lisandro Dalcin',
    'maintainer_email' : 'dalcinl@dalcinl.com',
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
PF_PKG_DIR   = abspath(join(PETSCFEM_DIR, '..', 'petscfem-packages'))

ANN      = PF_PKG_DIR
LIBRETTO = PF_PKG_DIR
MESCHACH = PF_PKG_DIR
NEWMAT   = PF_PKG_DIR
METIS    = PF_PKG_DIR
SSL      = PF_PKG_DIR
GLIB     = '/usr/include/glib-2.0'
GLIBCFG  = '/usr/lib/glib-2.0/include'


BOPT = env.get('BOPT') or 'g'
if PETSC_ARCH.endswith('O_c++') or PETSC_ARCH.endswith('O'):
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
                      join(PETSC_DIR, PETSC_ARCH, 'include'),
                      join(PETSC_DIR, 'bmake', PETSC_ARCH),
                      join(PETSC_DIR, 'include'),
                      GLIB, GLIBCFG,
                      join(ANN,       'include'),
                      join(LIBRETTO,  'include'),
                      join(NEWMAT,    'include/newmat'),
                      join(MESCHACH,  'include/meschach'),
                      join(METIS,     'include/metis'),
                      join(SSL,       'include/SSL'),
                      ],
    
    'libraries'    : ['petscksp', 'petscdm',
                      'petscmat', 'petscvec', 'petsc',
                      'mpich',
                      'ANN',  'ibretto', 'meschach',
                      'metis', 'newmat', 'simpleskts',
                      'glib-2.0', 'pthread',
                      ],
    
    'library_dirs' : [join(MPI_DIR,   'lib'),
                      join(PETSC_DIR, 'lib'),
                      join(PETSC_DIR, PETSC_ARCH, 'lib'),
                      join(PETSC_DIR, 'lib', PETSC_ARCH),
                      join(ANN,       'lib'),
                      join(LIBRETTO,  'lib'),
                      join(MESCHACH,  'lib'),
                      join(METIS,     'lib'),
                      join(NEWMAT,    'lib'),
                      join(SSL,       'lib'),
                      ],
    'runtime_library_dirs' : []
    }

def ext_modules(Extension):
    from os.path import abspath, join

    PETSCFEM_MACROS = config['define_macros']
    PETSCFEM_DIR = abspath('..')
    PETSCFEM_INCLUDE = [PETSCFEM_DIR, join(PETSCFEM_DIR, 'src')] \
                       + config['include_dirs']
    if BOPT[0] not in ('g', 'O'): raise SystemExit('invalid BOPT: %s' % BOPT)
    osfx = '_%s' % BOPT[0]
    PETSCFEM_LIBRARY = ['ns'       + osfx,
                        'advdif'   + osfx,
                        'petscfem' + osfx, ] * 2
    PETSCFEM_LIBRARY +=  config['libraries']
    PETSCFEM_LIBDIR = [join(PETSCFEM_DIR, 'src'),
                       join(PETSCFEM_DIR, 'applications/ns'),
                       join(PETSCFEM_DIR, 'applications/advdif'),
                       ] + config['library_dirs']
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
