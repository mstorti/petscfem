//__INSERT_LICENSE__
//$Id: spsolve.cpp,v 1.17.104.1 2007/02/19 20:23:56 mstorti Exp $

#include "sparse2.h"

//---:---<*>---:---<*>---:a---<*>---:---<*>---:---<*>---:---<*>---: 

#if PETSC_VERSION_(3,1,0) || PETSC_VERSION_(3,0,0)
inline PetscErrorCode VecDestroy_Compat(Vec *a)
{Vec b = *a; *a=0; return VecDestroy(b);}
#define VecDestroy VecDestroy_Compat
inline PetscErrorCode MatDestroy_Compat(Mat *a)
{Mat b = *a; *a=0; return MatDestroy(b);}
#define MatDestroy MatDestroy_Compat
inline PetscErrorCode KSPDestroy_Compat(KSP *a)
{KSP b = *a; *a=0; return KSPDestroy(b);}
#define KSPDestroy KSPDestroy_Compat
#endif

#if (PETSC_VERSION_MAJOR    == 2 && \
     PETSC_VERSION_MINOR    == 3 && \
     PETSC_VERSION_SUBMINOR == 3)
#if (PETSC_VERSION_RELEASE  == 1)
#define MatSetOption(mat,opt,flg) MatSetOption((mat),(opt))
#endif
#endif
//---:---<*>---:---<*>---:a---<*>---:---<*>---:---<*>---:---<*>---: 

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
    b = &*bb.begin();
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
#ifdef USE_SUPERLU
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

    assert(Bf==0);
    dCreate_Dense_Matrix(&B,m,1,b,m,DN,_D,GE); Bf=1;

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

    assert(Af==0);
    dCreate_CompCol_Matrix(&A,m,m,nnz,a,asub,xa,NR,_D,GE); Af=1;

    permc_spec = 2;
    get_perm_c(permc_spec, &A, perm_c);

    assert(Uf==0);
    assert(Lf==0);
    dgssv(&A, perm_c, perm_r, &L, &U, &B, &info); Uf=1; Lf=1;
    assert(info==0);

  }

      
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void SuperLUMat::solve_only() {

#if 0
    fact_and_solve();
    clean_factor();

#else
    int m,nnz,info;

    m = rows();
    assert(m == cols());
    nnz = size();

    assert(Bf==1);
    dCreate_Dense_Matrix(&B,m,1,b,m,DN,_D,GE); 

    assert(Uf=1);
    assert(Lf=1);
    dgstrs ("T",&L,&U,perm_r,perm_c,&B,&info); 
#endif
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

    assert(Af==1);
    assert(Bf==1);
    assert(Lf==1);
    assert(Uf==1);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    Af=0;
    Bf=0;
    Lf=0;
    Uf=0;
    
    delete[] perm_r;
    perm_r=NULL;
    delete[] perm_c;
    perm_c=NULL;
  }
#endif
    
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

    d_nnz_p = &*d_nnz.begin();

    e = end();
    for (j=0; j<m; j++) {
      row = find(j);
      assert(row != e);
      const Vec & v = row->second;
      d_nnz_p[j] = v.size();
    }

    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,m,m,
			   PETSC_NULL,d_nnz_p,&A); assert(!ierr);
    ierr =  MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
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

    // Create KSP and PETSc stuff
    ierr = KSPCreate(PETSC_COMM_SELF,&ksp); assert(!ierr); 
    ierr = KSPSetOperators(ksp,A,
			    A,SAME_NONZERO_PATTERN); assert(!ierr); 
    ierr = KSPGetPC(ksp,&pc); assert(!ierr); 
    
    ierr = KSPSetType(ksp,KSPPREONLY); assert(!ierr); 
    ierr = PCSetType(pc,PCLU); assert(!ierr); 
    ierr = PCFactorSetUseInPlace(pc); assert(!ierr);

    solve_only();

  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void PETScMat::solve_only() {

    int m,ierr,its;
    double *xx,*bb;
    ::Vec b_vec,x_vec;

    m = rows();
    assert(m == cols());

    ierr = VecCreateSeq(PETSC_COMM_SELF,m,&b_vec); assert(!ierr);
    ierr = VecDuplicate(b_vec,&x_vec); assert(!ierr);
    
    ierr = VecGetArray(b_vec,&bb); assert(!ierr); 
    memcpy(bb,b,m*sizeof(double));
    ierr = VecRestoreArray(b_vec,&bb); assert(!ierr); 
    
    ierr = KSPSolve(ksp,b_vec,x_vec); assert(!ierr); 
    ierr = KSPGetIterationNumber(ksp,&its); assert(!ierr);

    ierr = VecGetArray(x_vec,&xx); assert(!ierr); 
    memcpy(b,xx,m*sizeof(double));
    ierr = VecRestoreArray(x_vec,&xx); assert(!ierr); 

    ierr = VecDestroy(&x_vec); assert(!ierr); 
    ierr = VecDestroy(&b_vec); assert(!ierr); 

  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void PETScMat::clean_factor() {

    int ierr;
    ierr = MatDestroy(&A); assert(!ierr);
    ierr = KSPDestroy(&ksp); assert(!ierr);

  }
  
}

