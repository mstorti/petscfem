// -*- c++ -*-

%{
#include <petsc.h>
%}

%init {
if (!PetscInitializeCalled)
  PetscInitializeNoArguments();
}
