// -*- c++ -*-
// $Id: init.i,v 1.1.2.3 2006/03/28 22:13:25 rodrigop Exp $

%{
#include <petsc.h>
%}

%wrapper %{
static bool PetscFemBeganPetsc = false;
static void petscfem_atexit(void) {
  int flag;  MPI_Finalized(&flag); if (flag) return;
  if (!PetscFemBeganPetsc || PetscFinalizeCalled) return;
  PetscErrorCode ierr = PetscFinalize(); if (!ierr) return;
  fflush(stderr);
  fprintf(stderr, "PetscFinalize() failed [ierr: %d]\n",ierr);
  fflush(stderr);
}
%}

%init %{
  if (!PetscInitializeCalled && !PetscFinalizeCalled) {
    PetscInitializeNoArguments();
    PetscFemBeganPetsc = true;
    if (Py_AtExit(petscfem_atexit) < 0)
      PyErr_Warn(PyExc_RuntimeWarning, "cannot register "
		 "PetscFinalize() with Py_AtExit()");
  }
%}
