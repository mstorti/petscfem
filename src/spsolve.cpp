//__INSERT_LICENSE__
//$Id: spsolve.cpp,v 1.8 2001/11/08 03:28:12 mstorti Exp $

#include "sparse.h"

extern "C" {
  int PrintInt10(char *name, int len, int *x);
}

namespace Sparse {

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void Mat::solve(double *bb) {
    b = bb;
    fsm.solve();
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void Mat::solve(FullVec &bb) {
    assert(rows() == bb.length());
    b = bb.begin();
    fsm.solve();
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
  void SuperLUMat::fact_and_solve() {

    int j,m,nnz,curs,info;
    RowCIt row,e;
    VecCIt l,el;
    double *a;
    int *asub, *xa;
    int permc_spec;

    m = rows();
    assert(m == cols());
    nnz = size();

    dCreate_Dense_Matrix(&B,m,1,b,m,DN,_D,GE);

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
  void SuperLUMat::solve_only() {

    int j,m,nnz,curs,info;

    m = rows();
    assert(m == cols());
    nnz = size();

    dCreate_Dense_Matrix(&B,m,1,b,m,DN,_D,GE);

    dgstrs ("T",&L,&U,perm_r,perm_c,&B,&info);
      
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void SuperLUMat::clean_factor() {

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
    
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void SuperLUMat::fact_and_solve() {

    vector<int> d_nnz;
    int* d_nnz_p;
    RowCIt row,e;
    VecCIt l,el;

    m = rows();
    assert(m == cols());
    d_nnz.resize(m);

    d_nnz_p = d_nnz.begin();
    for (j=0; j<m; j++) {
      row = find(j);
      assert(row != e);
      const Vec & v = row->second;
      d_nnz_p[j] = v.size();
    }

    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,m,m,
			   PETSC_NULL,d_nnz_p,&A); assert(!ierr);
    ierr =  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR);
    for (j=0; j<m; j++) {
      row = find(j);
      assert(row != e);
      const Vec & v = row->second;
      el = v.end();
      for (l=v.begin(); l!=el; l++) {
	ierr = MatSetValues(A,1,&j,1,&l->first,
			    &l->second,INSERT_VALUES); assert(!ierr);
      }
    }
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); assert(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); assert(ierr);

  }

}

