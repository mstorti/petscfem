#!/bin/env python

# $id$

name    = 'petscfem4py'

version = '0.1.0'

classifiers = """
License :: GPL
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
finite elemets method
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

MPI_DIR      = env.get('MPI_DIR')    or '/usr/local/mpich2-1.0.2'
PETSC_DIR    = env.get('PETSC_DIR')  or '/usr/local/petsc/2.3.1'
PETSC_ARCH   = env.get('PETSC_ARCH') or 'linux-gnu-g_c++'
PETSCFEM_DIR = abspath('..')

SOFT_DIR     =  expanduser('~/PETSC')

GLIB1 = '/usr/include/glib-1.2'
GLIB2 = '/usr/lib/glib/include'

NEWMAT       = join(SOFT_DIR, 'NEWMAT/src')
LIBRETTO     = join(SOFT_DIR, 'libretto-2.1')
MESCHACH     = join(SOFT_DIR, 'meschach-1.2')
METIS        = join(SOFT_DIR, 'metis-4.0')
ANN_INCLUDE  = join(SOFT_DIR, 'ann/include')
ANN_LIBDIR   = join(SOFT_DIR, 'ann/lib')



# --------------------------------------------------------------------
# Extension modules
# --------------------------------------------------------------------

config = {
    'define_macros': [ ('MPICH_SKIP_MPICXX', None)],
    'include_dirs' : [join(MPI_DIR, 'include'),
                      join(PETSC_DIR, 'bmake', PETSC_ARCH),
                      join(PETSC_DIR, 'include'),
                      GLIB1, GLIB2, NEWMAT, LIBRETTO, MESCHACH, METIS, ANN_INCLUDE],
    
    'libraries'    : ['glib',
                      'mpich',
                      'petsc', 'petscvec', 'petscmat', 'petscksp',
                      'newmat',
                      'ibretto',
                      'mes',
                      'metis',
                      'ANN',
                      'simpleskts'],
    
    'library_dirs' : [join(MPI_DIR, 'lib'),
                      join(PETSC_DIR, 'lib', PETSC_ARCH),
                      '/usr/local/lib',expanduser('~/lib'),
                      NEWMAT, LIBRETTO, MESCHACH, METIS, ANN_LIBDIR],
    'runtime_library_dirs' : []
    }

def ext_modules(Extension):
    from os.path import abspath, join

    PETSCFEM_MACROS = config['define_macros']
    PETSCFEM_DIR = abspath('..')
    PETSCFEM_INCLUDE_DIR = [PETSCFEM_DIR, join(PETSCFEM_DIR, 'src')] \
                           + config['include_dirs']

    if 'g_c++' in PETSC_ARCH:
        PETSCFEM_LIBRARY = ['petscfem_g', 'ns_g'] * 2 + config['libraries']
    elif 'O_c++' in PETSC_ARCH:
        PETSCFEM_LIBRARY = ['petscfem_O', 'ns_O'] * 2 + config['libraries']
    else:
        raise SystemExit('invalid PETSC_ARCH: %s' % PETSC_ARCH)

    PETSCFEM_LIBRARY = ['petscfem_g', 'ns_g'] * 2 + config['libraries']
    PETSCFEM_LIBRARY_DIR = [join(PETSCFEM_DIR, 'src'),
                            join(PETSCFEM_DIR, 'applications/ns') ] \
                            + config['library_dirs']
    PETSCFEM_RT_LIB_DIR = PETSCFEM_LIBRARY_DIR \
                          + config['runtime_library_dirs']
    
    from glob import glob
    sources = glob('petscfem/petscfem.i') + glob('petscfem/*.cpp')
    #sources = glob('petscfem/petscfem.i') + glob('petscfem/NavierStokes.cpp')
    if 'petscfem/petscfem_wrap.cpp' in sources: 
        sources.remove('petscfem/petscfem_wrap.cpp')
    depends = glob('petscfem/*.i') + glob('petscfem/*.h')

    petscfem = Extension('petscfem4py._petscfem',
                         sources=sources,
                         depends=depends,
                         define_macros        = PETSCFEM_MACROS,
                         include_dirs         = PETSCFEM_INCLUDE_DIR,
                         libraries            = PETSCFEM_LIBRARY,
                         library_dirs         = PETSCFEM_LIBRARY_DIR,
                         runtime_library_dirs = PETSCFEM_RT_LIB_DIR,
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
    sources = glob('petscfem/*.cpp')
    if 'petscfem/petscfem_wrap.cpp' in sources: 
        sources.remove('petscfem/petscfem_wrap.cpp')

    return [('petscfem', dict(sources      = sources,
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
            self.swigflags.append('-O')
            self.swigflags.append('-nofastproxy')
            self.swigflags.append('-noh')
            
    setup(packages     = ['petscfem4py'],
    	  package_dir  = {'petscfem4py' : 'petscfem'},
          ext_modules  = ext_modules(Extension),
          #libraries    = c_libraries(),
          cmdclass     = {'build_src' : build_src},
          **metadata)

    
if __name__ == '__main__':

    from distutils import sysconfig
    cvars = sysconfig.get_config_vars()
    cflags = cvars['OPT'].split()
    cflags.remove('-Wall')
    cvars['OPT'] = str.join(' ', cflags) 
    setup()

# --------------------------------------------------------------------
