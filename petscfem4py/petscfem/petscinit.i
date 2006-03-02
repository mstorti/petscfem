// -*- c++ -*-
// $Id: petscinit.i,v 1.1.2.1 2006/03/02 21:37:12 rodrigop Exp $

%{
#include <petsc.h>
%}

%wrapper %{
static bool PetscFemBeganPetsc = false;
static void petscfem_atexit(void) {
  if (!PetscFemBeganPetsc || PetscFinalizeCalled) return;
  PetscErrorCode ierr = PetscFinalize(); if (!ierr) return;
  fflush(stderr);
  fprintf(stderr, "PetscFinalize() failed [ierr: %d]\n",ierr);
  fflush(stderr);
}
%}

%init {
  if (!PetscInitializeCalled && !PetscFinalizeCalled) {
    PetscInitializeNoArguments();
    PetscFemBeganPetsc = true;
    if (Py_AtExit(petscfem_atexit) < 0)
      PyErr_Warn(PyExc_RuntimeWarning, "cannot register "
		 "PetscFinalize() with Py_AtExit()");
  }
}
