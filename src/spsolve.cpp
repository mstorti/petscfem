//__INSERT_LICENSE__
//$Id: spsolve.cpp,v 1.6 2001/09/29 00:57:38 mstorti Exp $

#include "sparse.h"

extern "C" {
  int PrintInt10(char *name, int len, int *x);
}

namespace Sparse {

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void Mat::solve(FullVec &bb) {
    b = &bb;
    fsm.solve();
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void Mat::fact_and_solve() {

    int j,m,nnz,curs,info;
    RowCIt row,e;
    VecCIt l,el;
    double *a;
    int *asub, *xa;
    int permc_spec;

    m = rows();
    assert(m == cols());
    assert(m == b->length());
    nnz = size();

    dCreate_Dense_Matrix(&B,m,1,b->begin(),m,DN,_D,GE);

    a = new double[nnz];
    asub = new int[nnz];
    xa = new int[m+1];
    perm_r = new int[m];
    perm_c = new int[m];
      
    xa[0] = 0;
    e = end();
    curs = 0;
    for (j=0; j<m; j++) {
      row = find(j);
      assert(row != e);
      const Vec & v = row->second;
      xa[j+1] = xa[j] + v.size();
      el = v.end();
      for (l=v.begin(); l!=el; l++) {
	a[curs] = l->second;
	asub[curs] = l->first;
	curs++;
      }
    }

    dCreate_CompCol_Matrix(&A,m,m,nnz,a,asub,xa,NR,_D,GE);

    permc_spec = 2;
    get_perm_c(permc_spec, &A, perm_c);

    dgssv(&A, perm_c, perm_r, &L, &U, &B, &info);
    assert(info==0);

  }

      
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void Mat::solve_only() {

    int j,m,nnz,curs,info;

    m = rows();
    assert(m == cols());
    assert(m == b->length());
    nnz = size();

    dCreate_Dense_Matrix(&B,m,1,b->begin(),m,DN,_D,GE);

    dgstrs ("T",&L,&U,perm_r,perm_c,&B,&info);
      
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  Mat & Mat::clear() {

    fsm.clear();
    return *this;

  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void Mat::clean_mat() {
    
    map<int,Vec>::clear();

  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void Mat::clean_factor() {

#if 0
    dPrint_CompCol_Matrix("A", &A);
    dPrint_CompCol_Matrix("U", &U);
    dPrint_SuperNode_Matrix("L", &L);
    PrintInt10("\nperm_r", m, perm_r);
    delete[] asub;
    delete[] xa;
    delete[] a;
#endif

    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    
    delete[] perm_r;
    delete[] perm_c;
  }
    
}

