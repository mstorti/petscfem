## -*- mode: cmake -*-
cmake_minimum_required(VERSION 3.15)
project(PETScFEM VERSION 1.0.0 LANGUAGES CXX Fortran C)

# Configuración de C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/../petscfem-opt)

set(BUILD_SHARED_LIBS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Esto desactiva el "New DTags" (RUNPATH) y vuelve al RPATH clásico
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--disable-new-dtags")

find_package(PkgConfig REQUIRED)

add_compile_definitions(USE_DLEF)
add_compile_definitions(USE_HDF5)

cmake_host_system_information(RESULT HOSTNAME QUERY HOSTNAME)
if(${HOSTNAME} STREQUAL galileo)
  set(PFPCKGDIR $ENV{HOME}/PETSC/petscfem-packages)
  set(NEWMAT_ROOT ${PFPCKGDIR}) 
  set(LIBRETTO_ROOT ${PFPCKGDIR}) 
  set(METIS_ROOT ${PFPCKGDIR}) 
  set(MESCHACH_ROOT ${PFPCKGDIR}) 
  set(HDF5_ROOT $ENV{HOME}/SOFT/hdf5-1814)
  set(PETSC_DIR "/home/mstorti/PETSC/petsc-3.2-p7")
  set(PETSC_ARCH "linux-gnu-O")
elseif(${HOSTNAME} STREQUAL urubu)
  set(PUB /share/storage/pub/mstorti)
  set(PFPCKGDIR ${PUB}/petscfem-packages)
  set(NEWMAT_ROOT ${PFPCKGDIR}) 
  set(LIBRETTO_ROOT ${PFPCKGDIR}) 
  set(METIS_ROOT ${PFPCKGDIR}) 
  set(MESCHACH_ROOT ${PFPCKGDIR}) 
  set(HDF5_ROOT ${PUB}/hdf5-1814)
  set(PETSC_DIR ${PUB}/petsc-3.2-p7)
  set(PETSC_ARCH linux-gnu-O)
elseif(${HOSTNAME} STREQUAL pirayu)
  set(PUB /home/pub/mstorti)
  set(PFPCKGDIR ${PUB}/petscfem-packages)
  set(NEWMAT_ROOT ${PFPCKGDIR}) 
  set(LIBRETTO_ROOT ${PFPCKGDIR}) 
  set(METIS_ROOT ${PFPCKGDIR}) 
  set(MESCHACH_ROOT ${PFPCKGDIR}) 
  set(HDF5_ROOT ${PUB}/hdf5-1814)
  set(PETSC_DIR ${PUB}/petsc-3.2-p7)
  set(PETSC_ARCH linux-gnu-O)
  list(APPEND CMAKE_PREFIX_PATH ${PUB}/jsoncpp)
  set(jsoncpp_DIR /home/pub/mstorti/jsoncpp/lib64/cmake/jsoncpp)
#  find_package(jsoncpp REQUIRED CONFIG)
endif()

set(PETSC_ROOT ${PETSC_DIR}/${PETSC_ARCH})
add_compile_definitions(PETSC_BOPT_FLAG="${PETSC_ARCH}")

## HDF5
add_library(hdf5_lib INTERFACE)
list(APPEND CMAKE_INSTALL_RPATH "/home/mstorti/SOFT/hdf5-1814/lib")
link_directories(${HDF5_ROOT}/lib)
target_include_directories(hdf5_lib INTERFACE ${HDF5_ROOT}/include)
find_library(H5_STUB_LIB NAMES hdf5_f90cstub
  PATHS ${HDF5_ROOT}/lib NO_DEFAULT_PATH REQUIRED)
target_link_libraries(hdf5_lib INTERFACE hdf5 hdf5_cpp ${H5_STUB_LIB})

## NEWMAT
add_library(newmat_lib INTERFACE)
target_include_directories(newmat_lib INTERFACE ${NEWMAT_ROOT}/include/newmat)
find_library(NEWMAT_BINARY NAMES newmat PATHS ${NEWMAT_ROOT}/lib REQUIRED)
target_link_libraries(newmat_lib INTERFACE ${NEWMAT_BINARY})

## LIBRETTO
add_library(libretto_lib INTERFACE)
target_include_directories(libretto_lib INTERFACE ${LIBRETTO_ROOT}/include)
find_library(LIBRETTO_BINARY NAMES ibretto PATHS ${LIBRETTO_ROOT}/lib REQUIRED)
target_link_libraries(libretto_lib INTERFACE ${LIBRETTO_BINARY})

## METIS
add_library(metis_lib INTERFACE)
target_include_directories(metis_lib INTERFACE ${METIS_ROOT}/include/metis)
find_library(METIS_BINARY NAMES metis PATHS ${METIS_ROOT}/lib
  NO_DEFAULT_PATH REQUIRED)
target_link_libraries(metis_lib INTERFACE ${METIS_BINARY})
list(APPEND CMAKE_INSTALL_RPATH "${METIS_ROOT}/lib")

## MESCHACH
add_library(meschach_lib INTERFACE)
target_include_directories(meschach_lib INTERFACE ${MESCHACH_ROOT}/include/meschach)
find_library(MESCHACH_BINARY NAMES meschach PATHS ${MESCHACH_ROOT}/lib REQUIRED)
target_link_libraries(meschach_lib INTERFACE ${MESCHACH_BINARY})

## GLIB
pkg_check_modules(GLIB REQUIRED glib-2.0)

## MPI
find_package(MPI)

## X (for PETSc)
find_package(X11 REQUIRED)

## PETSC
find_library(PETSC_LIBRARY_FILE
    NAMES petsc
    PATHS ${PETSC_ROOT}/lib
    NO_DEFAULT_PATH
    REQUIRED)
add_library(petsc_lib INTERFACE)
target_include_directories(petsc_lib INTERFACE
  ${PETSC_DIR}/include
  ${PETSC_ROOT}/include)
target_link_libraries(petsc_lib INTERFACE
  ${PETSC_LIBRARY_FILE}
  X11::X11
  MPI::MPI_C
  MPI::MPI_Fortran
  ${CMAKE_DL_LIBS})

## BLAS LAPACK
pkg_check_modules(FLEXIBLAS QUIET flexiblas)
add_library(blas_lapack_lib INTERFACE)

if(FLEXIBLAS_FOUND)
  message(STATUS "FlexiBLAS detectado (Fedora moderno)")
  target_link_libraries(blas_lapack_lib INTERFACE ${FLEXIBLAS_LIBRARIES})
  target_include_directories(blas_lapack_lib INTERFACE ${FLEXIBLAS_INCLUDE_DIRS})
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  target_link_libraries(blas_lapack_lib INTERFACE
    ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  message(STATUS "--- DEBUG BLAS/LAPACK ---")
  message(STATUS "BLAS Libs: ${BLAS_LIBRARIES}")
  message(STATUS "LAPACK Libs: ${LAPACK_LIBRARIES}")
endif()

## FLEXIBLAS
# pkg_check_modules(FLEXIBLAS REQUIRED flexiblas)
# add_library(blas_lib INTERFACE)
# target_include_directories(blas_lib INTERFACE ${FLEXIBLAS_INCLUDE_DIRS})
# target_link_libraries(blas_lib INTERFACE ${FLEXIBLAS_LIBRARIES})
# target_compile_options(blas_lib INTERFACE ${FLEXIBLAS_CFLAGS_OTHER})

## BLAS
# add_library(blas_lib INTERFACE)
# target_include_directories(blas_lib INTERFACE /usr/include/flexiblas)
# link_directories(/usr/lib64)
# target_link_libraries(blas_lib INTERFACE flexiblas lapacke)

## JSONCPP
find_package(jsoncpp REQUIRED CONFIG)
## pkg_check_modules(JSONCPP REQUIRED jsoncpp)
## add_library(jsoncpp_lib INTERFACE)
target_include_directories(jsoncpp_lib INTERFACE ${JSONCPP_INCLUDE_DIRS})
target_link_libraries(jsoncpp_lib INTERFACE ${JSONCPP_LIBRARIES})
# (Opcional) Si necesitas definir flags específicos
## target_compile_options(jsoncpp_lib INTERFACE ${JSONCPP_CFLAGS_OTHER})

# Incluir subdirectorios
add_subdirectory(src)
add_subdirectory(applications/advdif)
add_subdirectory(applications/ns)
## add_subdirectory(applications/ns)

# Agregamos todas las librerías auxiliares al conjunto de exportación "petscfemTargets"
install(TARGETS
  petsc_lib
  hdf5_lib 
  newmat_lib 
  libretto_lib 
  metis_lib 
  meschach_lib
  blas_lapack_lib
  ## jsoncpp_lib
  EXPORT petscfemTargets)
