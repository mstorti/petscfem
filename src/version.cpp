// -*- mode: c++ -*-
#include <petsc.h>

void print_copyright() { 
  PetscPrintf(PETSC_COMM_WORLD,
"========================================================================
This is PETSc-FEM, version 2.28 (Mon Dec 25 21:59:34 2000 localhost.localdomain)
Copyright (C) 1999,2000 Mario A. Storti.  
PETSc-FEM comes with ABSOLUTELY NO WARRANTY.  This is free software,
and you are welcome to redistribute it under certain conditions.
========================================================================

");
}
