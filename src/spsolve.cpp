//__INSERT_LICENSE__
//$Id: spsolve.cpp,v 1.13 2001/11/13 17:34:25 mstorti Exp $

#include "sparse2.h"

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
  void PETScMat::fact_and_solve() {

    vector<int> d_nnz;
    int *d_nnz_p,m,j,k,ierr;
    RowCIt row,e;
    VecCIt l,el;
    double w;

    m = rows();
    assert(m == cols());
    d_nnz.resize(m);

    d_nnz_p = d_nnz.begin();

    e = end();
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
	k = l->first;
	w = l->second;
	ierr = MatSetValues(A,1,&j,1,&k,&w,INSERT_VALUES); assert(!ierr);
      }
    }
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); assert(!ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); assert(!ierr);

#if 0

    RowCIt row,e;
    VecCIt l,el;

    m = rows();
    assert(m == cols());
    nnz = size();
#endif

    // Create SLES and PETSc stuff
    ierr = SLESCreate(PETSC_COMM_SELF,&sles); assert(!ierr); 
    ierr = SLESSetOperators(sles,A,
			    A,SAME_NONZERO_PATTERN); assert(!ierr); 
    ierr = SLESGetKSP(sles,&ksp); assert(!ierr); 
    ierr = SLESGetPC(sles,&pc); assert(!ierr); 
    
    ierr = KSPSetType(ksp,KSPPREONLY); assert(!ierr); 
    ierr = PCSetType(pc,PCLU); assert(!ierr); 
    ierr = PCLUSetUseInPlace(pc); assert(!ierr);

    solve_only();

  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void PETScMat::solve_only() {

    int m,j,ierr,k,its;
    double *xx,*bb;
    ::Vec b_vec,x_vec;

    m = rows();
    assert(m == cols());

    ierr = VecCreateSeq(PETSC_COMM_SELF,m,&b_vec); assert(!ierr);
    ierr = VecDuplicate(b_vec,&x_vec); assert(!ierr);
    
    ierr = VecGetArray(b_vec,&bb); assert(!ierr); 
    memcpy(bb,b,m*sizeof(double));
    ierr = VecRestoreArray(b_vec,&bb); assert(!ierr); 
    
    ierr = SLESSolve(sles,b_vec,x_vec,&its); assert(!ierr); 

    ierr = VecGetArray(x_vec,&xx); assert(!ierr); 
    memcpy(b,xx,m*sizeof(double));
    ierr = VecRestoreArray(x_vec,&xx); assert(!ierr); 

    ierr = VecDestroy(x_vec); assert(!ierr); 
    ierr = VecDestroy(b_vec); assert(!ierr); 

  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void PETScMat::clean_factor() {

    int ierr;
    ierr = MatDestroy(A); assert(!ierr);
    ierr = SLESDestroy(sles); assert(!ierr);

  }
  
}

