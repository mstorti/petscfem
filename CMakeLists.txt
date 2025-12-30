## -*- mode: cmake -*-
cmake_minimum_required(VERSION 3.15)
project(PETScFEM LANGUAGES CXX)

# Configuración de C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(PFPCKGDIR $ENV{HOME}/PETSC/petscfem-packages)
set(NEWMAT_ROOT ${PFPCKGDIR}) 

# 2. Crear un objetivo (target) para Newmat
add_library(newmat_lib INTERFACE)

# 3. Definir dónde están los headers y la librería (.a o .so)
target_include_directories(newmat_lib INTERFACE ${NEWMAT_ROOT}/include/newmat)

# Buscamos el archivo físico de la librería para no hardcodear nombres
find_library(NEWMAT_BINARY NAMES newmat PATHS ${NEWMAT_ROOT}/lib REQUIRED)
target_link_libraries(newmat_lib INTERFACE ${NEWMAT_BINARY})

# Buscar dependencias (PETSc es clave aquí)
find_package(PkgConfig REQUIRED)
## pkg_check_modules(PETSC REQUIRED petsc)

set(PFPCKGDIR $ENV{HOME}/PETSC/petscfem-packages)
set(HDF5_ROOT /home/mstorti/SOFT/hdf5-1814)
set(CS_ROOT /home/scratch/mstorti/TO-BACKUP/DATA-CD/saturne-opt)
set(PETSC_DIR /home/mstorti/PETSC/petsc-3.2-p7)
set(PETSC_ARCH linux-gnu-O)
set(PETSC_ROOT ${PETSC_DIR}/${PETSC_ARCH})
set(PETSC_INCLUDES ${PETSC_DIR}/include ${PETSC_DIR}/${PETSC_ARCH}/include)

# Incluir subdirectorios
add_subdirectory(src)
## add_subdirectory(applications/advdif)
## add_subdirectory(applications/ns)

