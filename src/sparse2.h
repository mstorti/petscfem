// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: sparse2.h,v 1.3 2002/09/05 18:23:52 mstorti Exp $
#ifndef SPARSE2_H
#define SPARSE2_H

#include <src/sparse.h>
#include <petscsles.h>

using namespace Random;

namespace Sparse {

  class SuperLUMat : public Mat {
  private:
    /// Factored matrix
    SuperMatrix A,L,U,B;
    /// Flags that indicate if the matrix was allocated or not
    int Af, Lf, Uf, Bf;
    int *perm_r, *perm_c, *etree;
  public:
    SuperLUMat() : Af(0), Lf(0), Uf(0), Bf(0), 
      perm_r(NULL), perm_c(NULL) {}
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
