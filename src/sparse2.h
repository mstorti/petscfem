// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: sparse2.h,v 1.1 2001/11/09 03:05:42 mstorti Exp $
#ifndef SPARSE2_H
#define SPARSE2_H

#include <src/sparse.h>
#include <sles.h>

using namespace Random;

namespace Sparse {

  class SuperLUMat : public Mat {
  private:
    /// Factored matrix
    SuperMatrix A,L,U,B;
    int *perm_r, *perm_c;
  public:
    void clean_factor();
    void solve_only();
    void fact_and_solve();
  };

  class PETScMat : public Mat {
  private:
    /// Factored matrix
    ::Mat A;
    SLES sles;
    KSP ksp;
    PC pc;
  public:
    void clean_factor();
    void solve_only();
    void fact_and_solve();
  };

}

#endif
