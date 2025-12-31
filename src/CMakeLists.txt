## -*- mode: cmake -*-

# src/CMakeLists.txt
## file(GLOB_RECURSE SRC_FILES *.cpp)
add_library(petscfem 
  amplidl.cpp arglist.cpp autostr.cpp blessf.cpp cloud.cpp cloud2.cpp
  debug.cpp distmat.cpp dlhook.cpp dofmap.cpp dofmap2.cpp dofmap3.cpp
  dvecpar.cpp dvector.cpp dxelmst.cpp dxhook.cpp elemlist.cpp
  elemset.cpp elmsupl.cpp elmsupl2.cpp fastmat.cpp fastmat2.cpp
  fastmat2b.cpp fastmat2c.cpp fastmat2d.cpp fem.cpp fm2ctx2.cpp
  fm2eperl.cpp fm2eperl2.cpp fm2temp.cpp fstack.cpp gatherer.cpp
  generror.cpp getarray.cpp getprop.cpp gpdata.cpp graph.cpp
  graphdv.cpp hook.cpp idmap.cpp iisdcr.cpp iisdgraph.cpp iisdmat.cpp
  iisdmat2.cpp interpola.cpp lagmul.cpp linkgraph.cpp mainutl.cpp
  matchg.cpp measperf.cpp metisprt.cpp penalize.cpp petscmat.cpp
  pfmat.cpp pfobject.cpp propef.cpp randomg.cpp readmesh.cpp
  rhhook.cpp secant.cpp shllhook.cpp sparse.cpp spdirect.cpp
  spetscmat.cpp spsolve.cpp srfgath.cpp state.cpp sttfilter.cpp
  surf2vol.cpp syncbuff.cpp tempfun.cpp tempfun2.cpp texthash.cpp
  texthf.cpp timestat.cpp util2.cpp util3.cpp utils.cpp fm2prod2.cpp
  fmgemmcode.cpp h5utils.cpp version.cpp dgefa.f seval.f spline.f
  getarrgr.tab.c)

target_include_directories(petscfem PUBLIC ${CMAKE_SOURCE_DIR})

# Vincular con PETSc
target_compile_options(petscfem PUBLIC ${PETSC_CFLAGS_OTHER})

target_link_libraries(petscfem PUBLIC
  newmat_lib
  libretto_lib
  metis_lib
  meschach_lib
  MPI::MPI_CXX
  ${PETSC_LIBRARIES}
  ${GLIB_LIBRARIES}
  hdf5_lib
  blas_lib
  jsoncpp_lib )

# IMPORTANTE: Debes agregar los directorios de include explícitamente si usas variables
target_include_directories(petscfem PUBLIC ${GLIB_INCLUDE_DIRS})
