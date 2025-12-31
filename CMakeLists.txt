## -*- mode: cmake -*-
cmake_minimum_required(VERSION 3.15)
project(PETScFEM LANGUAGES CXX)

# Configuración de C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(PFPCKGDIR $ENV{HOME}/PETSC/petscfem-packages)
set(NEWMAT_ROOT ${PFPCKGDIR}) 

set(PETSC_DIR /home/mstorti/PETSC/petsc-3.2-p7)
set(PETSC_ARCH linux-gnu-O)

set(PETSC_ROOT ${PETSC_DIR}/${PETSC_ARCH})
set(PETSCLIB -L${PETSC_ROOT}/lib petsc)

# 2. Crear un objetivo (target) para Newmat
add_library(newmat_lib INTERFACE)
target_include_directories(newmat_lib INTERFACE ${NEWMAT_ROOT}/include/newmat ${PETSC_INCLUDES})
find_library(NEWMAT_BINARY NAMES newmat PATHS ${NEWMAT_ROOT}/lib REQUIRED)
target_link_libraries(newmat_lib INTERFACE ${NEWMAT_BINARY})

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

