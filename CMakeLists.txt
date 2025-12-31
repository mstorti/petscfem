## -*- mode: cmake -*-
cmake_minimum_required(VERSION 3.15)
project(PETScFEM LANGUAGES CXX)

# Configuración de C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_definitions(USE_DLEF)
add_compile_definitions(USE_HDF5)

set(PFPCKGDIR $ENV{HOME}/PETSC/petscfem-packages)
set(NEWMAT_ROOT ${PFPCKGDIR}) 
set(LIBRETTO_ROOT ${PFPCKGDIR}) 
set(METIS_ROOT ${PFPCKGDIR}) 
set(MESCHACH_ROOT ${PFPCKGDIR}) 
set(HDF5ROOT ${HOME}/SOFT/hdf5-1814)

set(PETSC_DIR /home/mstorti/PETSC/petsc-3.2-p7)
set(PETSC_ARCH linux-gnu-O)

set(PETSC_ROOT ${PETSC_DIR}/${PETSC_ARCH})

add_library(hdf5_lib INTERFACE)
target_include_directories(hdf5_lib INTERFACE ${HDF5_ROOT}/include/hdf5)
target_link_libraries(hdf5_lib INTERFACE hdf5 hdf5_cpp)

add_library(newmat_lib INTERFACE)
target_include_directories(newmat_lib INTERFACE ${NEWMAT_ROOT}/include/newmat)
find_library(NEWMAT_BINARY NAMES newmat PATHS ${NEWMAT_ROOT}/lib REQUIRED)
target_link_libraries(newmat_lib INTERFACE ${NEWMAT_BINARY})

add_library(libretto_lib INTERFACE)
target_include_directories(libretto_lib INTERFACE ${LIBRETTO_ROOT}/include)
find_library(LIBRETTO_BINARY NAMES ibretto PATHS ${LIBRETTO_ROOT}/lib REQUIRED)
target_link_libraries(libretto_lib INTERFACE ${LIBRETTO_BINARY})

add_library(metis_lib INTERFACE)
target_include_directories(metis_lib INTERFACE ${METIS_ROOT}/include/metis)
find_library(METIS_BINARY NAMES metis PATHS ${METIS_ROOT}/lib REQUIRED)
target_link_libraries(metis_lib INTERFACE ${METIS_BINARY})

# MESCHACH_INC_DIR = $(MESCHACH)/include/meschach
# MESCHACH_LIB_DIR = $(MESCHACH)/lib
# MESCHACH_LIBNAME = meschach

add_library(meschach_lib INTERFACE)
target_include_directories(meschach_lib INTERFACE ${MESCHACH_ROOT}/include/meschach)
find_library(MESCHACH_BINARY NAMES meschach PATHS ${MESCHACH_ROOT}/lib REQUIRED)
target_link_libraries(meschach_lib INTERFACE ${MESCHACH_BINARY})

find_package(PkgConfig REQUIRED)
pkg_check_modules(GLIB REQUIRED glib-2.0)

# Esto define automáticamente las variables:
# ${GLIB_INCLUDE_DIRS} -> Tiene ambas rutas (/usr/include... y /usr/lib64/...)
# ${GLIB_LIBRARIES}    -> Tiene la librería (glib-2.0)

find_package(MPI)

link_directories(${PETSC_ROOT}/lib)
link_libraries(petsc)
set(PETSC_INCLUDES ${PETSC_DIR}/include ${PETSC_DIR}/${PETSC_ARCH}/include)
include_directories(${PETSC_INCLUDES})

set(PFPCKGDIR $ENV{HOME}/PETSC/petscfem-packages)
set(HDF5_ROOT /home/mstorti/SOFT/hdf5-1814)
set(CS_ROOT /home/scratch/mstorti/TO-BACKUP/DATA-CD/saturne-opt)

# Incluir subdirectorios
add_subdirectory(src)
## add_subdirectory(applications/advdif)
## add_subdirectory(applications/ns)

