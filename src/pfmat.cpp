//__INSERT_LICENSE__
//$Id: pfmat.cpp,v 1.1.2.1 2002/01/12 15:19:00 mstorti Exp $

#include <mat.h>

#include <src/fem.h>
#include <src/utils.h>
#include <src/dofmap.h>
#include <src/elemset.h>
// #pragma implementation "PFMat"
#include <src/pfmat.h>
#include <src/pfptscmat.h>
#include <src/iisdmat.h>
#include <src/petscmat.h>
#include <src/spdirect.h>

PFMat * PFMat::dispatch(int N,DofPartitioner &part,const char *s) {
  PFMat *A;
  IISDMat *AA;
  // IISD solver with PETSc or SuperLU local solver
  if (!strcmp(s,"iisd_superlu")) {
    AA =  new IISDMat(N,N,part,PETSC_COMM_WORLD);
    AA->local_solver = IISDMat::SuperLU;
    return AA;
  } else if (!strcmp(s,"iisd_petsc")) {
    AA =  new IISDMat(N,N,part,PETSC_COMM_WORLD);
    AA->local_solver = IISDMat::PETSc;
    return AA;
  } else if (!strcmp(s,"iisd")) {
    // local solver is chosen by default
    AA =  new IISDMat(N,N,part,PETSC_COMM_WORLD);
    return AA;
  } else if (!strcmp(s,"petsc")) {
    // PETSc (iterative) solver 
    A = new PETScMat(N,N,part,PETSC_COMM_WORLD);
    return A;
  } else if (!strcmp(s,"direct_superlu")) {
    A = new SparseDirect(N,"SuperLU");
    return A;
  } else if (!strcmp(s,"direct") || 
	     !strcmp(s,"direct_petsc")) {
    A = new SparseDirect(N,"PETSc");
    return A;
  } else {
    PETSCFEM_ERROR("PFMat type not known: %s\n",s);
  }
}

