// -*- c++ -*-
// $Id$

%{
#include <petsc.h>
%}

%header %{
#if defined(PETSC_HAVE_MPIUNI)
#if !defined(MPI_Finalized)
#define MPI_Finalized(flag) (((flag)?(*(flag)=0):0),0)
#endif
#endif
%}

%wrapper %{
static bool PetscFemBeganPetsc = false;
static void petscfem_atexit(void) {
  int flag;  MPI_Finalized(&flag); if (flag) return;
  if (!PetscFemBeganPetsc || PetscFinalizeCalled) return;
  PetscErrorCode ierr = PetscFinalize(); 
  if (ierr) {
    fflush(stderr);
    fprintf(stderr, "PetscFinalize() failed [ierr: %d]\n",ierr);
    fflush(stderr);
  }
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
