// -*- mode: c++ -*-
#include <petsc.h>

void print_copyright() { 
  PetscPrintf(PETSC_COMM_WORLD,
"========================================================================
This is PETSc-FEM, version <:=`cat ../VERSION`:>
Copyright (C) 1999,2000 Mario A. Storti.  
PETSc-FEM comes with ABSOLUTELY NO WARRANTY.  This is free software,
and you are welcome to redistribute it under certain conditions.
========================================================================

");
}
