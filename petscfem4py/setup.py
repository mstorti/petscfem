#!/bin/env python

metadata = {'name' : 'petscfem4py',
            'version': '0.1.0',
            }

# --------------------------------------------------------------------
# Extension modules
# --------------------------------------------------------------------

config = {
    'define_macros': [ ('MPICH_SKIP_MPICXX', None)],
    'include_dirs' : ['/usr/include/glib-1.2', '/usr/lib/glib/include',
                      '/usr/local/mpich2-1.0.2/include',
                      '/u/rodrigop/PETSC/petsc-2.3.1-p0',
                      '/u/rodrigop/PETSC/petsc-2.3.1-p0/bmake/linux-gnu-g_c++',
                      '/u/rodrigop/PETSC/petsc-2.3.1-p0/include',
                      '/u/rodrigop/PETSC/NEWMAT/src',
                      '/u/rodrigop/PETSC/libretto-2.1',
                      '/u/rodrigop/PETSC/meschach-1.2',
                      '/u/rodrigop/PETSC/metis-4.0'
                      ],
    
    'library_dirs' : ['/usr/local/mpich2-1.0.2/lib',
                      '/u/rodrigop/PETSC/petsc-2.3.1-p0/lib/linux-gnu-g_c++',
                      '/u/rodrigop/PETSC/NEWMAT/src',
                      '/usr/local/lib',
                      '/u/rodrigop/PETSC/meschach-1.2',
                      '/u/rodrigop/PETSC/metis-4.0',
                      '/u/rodrigop/PETSC/ann/lib',
                      '/u/rodrigop/lib',
                      ],
    'libraries'    : ['glib',
                      'mpich',
                      'petsc', 'petscvec', 'petscmat', 'petscksp',
                      'newmat',
                      'ibretto',
                      'mes',
                      'metis',
                      'ANN',
                      'simpleskts'],
    
    }

def ext_modules(Extension):
    from os.path import abspath, join

    PETSCFEM_MACROS = config['define_macros']
    PETSCFEM_DIR = abspath('..')
    PETSCFEM_INCLUDE_DIR = [ join(PETSCFEM_DIR, d) for d in 
                            ['', 'src',
			     ]
			   ] + config['include_dirs']
    
    PETSCFEM_LIBRARY = ['ns_g', 'petscfem_g', 'ns_g', 'petscfem_g'] + config['libraries']
    PETSCFEM_LIBRARY_DIR = [join(PETSCFEM_DIR, 'src'),
                            join(PETSCFEM_DIR, 'applications/ns'),
                            ] + config['library_dirs']
    from glob import glob
    sources = glob('petscfem/petscfem.i') + glob('petscfem/*.cpp')
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
                         runtime_library_dirs = PETSCFEM_LIBRARY_DIR,
                         language='c++',
                       )
    return [petscfem]

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
            
    setup(packages     = ['petscfem4py'],
    	  package_dir  = {'petscfem4py' : 'petscfem'},
          ext_modules  = ext_modules(Extension),
          cmdclass     = {'build_src' : build_src},
          **metadata)

    
if __name__ == '__main__':
    setup()

# --------------------------------------------------------------------
