// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: monitor2.h,v 1.1 2003/07/07 21:15:26 mstorti Exp $
#ifndef PETSCFEM_MONITOR2_H
#define PETSCFEM_MONITOR2_H

#include <cstdio>
#include <petsc.h>
#include <src/monitor.h>

class DefaultMonitor : public Monitor {
 private:
  PFPETScMat *A;
  friend class PFPETScMat;
 public:
  void step(int n,double rnorm) { 
    if (A->print_internal_loop_conv) {
      if (n==0) PetscPrintf(A->comm,
			    " Begin internal iterations "
			    "--------------------------------------\n");
      PetscPrintf(A->comm,
		  "iteration %d KSP Residual_norm = %14.12e \n",n,rnorm);
    }
  }
};

#endif
