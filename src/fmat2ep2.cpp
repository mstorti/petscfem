//<require "prep.pl"; require "fmat2.pl"; //>//
//<=$warn_dont_modify //>

//__INSERT_LICENSE__
//$Id: fmat2ep2.cpp,v 1.1 2002/11/28 15:13:25 mstorti Exp $
#include <math.h>
#include <stdio.h>

#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fastlib2.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class eig_cache : public FastMatSubCache {
public:
  const double **A;
  double **W;
  double *A_c,*work,*W_c;
  int lwork;
  eig_cache() { W=NULL; A=NULL; A_c=W_c=NULL; work=NULL; }
  ~eig_cache() { 
    delete[] A;
    delete[] W;
    delete[] A_c;
    delete[] W_c;
    delete[] work;
  }
};

FastMat2 & FastMat2::eig(const FastMat2 & A) { 
  assert(0); // Not implemented yet
  return *this; 
}

FastMat2 & 
FastMat2::eig(const FastMat2 & A, FastMat2 &V,int compute_eigen_vectors=1) { 
  assert(0); // Not implemented yet
  return *this; 
}

FastMat2 & 
FastMat2::seig(const FastMat2 & A) {
  seig(A,(FastMat2 & )A,0);
  return *this; 
}

extern "C" void dsyev_(const char*jobz,const char *uplo,int *n,double *a,int *lda,
		       double *w,double *work, int *lwork, int *info);

//<$seig=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 &
FastMat2::seig(const FastMat2 & A, FastMat2 &V,int compute_eigen_vectors=1) { 

  __CACHE_OPERATIONS__;

  eig_cache *ecache;
  if (!was_cached) {
    assert(compute_eigen_vectors==0);
    Indx Adims;
    A.get_dims(Adims);
    int ndims = Adims.size();
    assert (ndims==2);
    assert (Adims[0] == Adims[1]);
    int m= Adims[0];

    Indx Wdims,Wdims_c;
    Wdims_c.push_back(m);
    if (!defined) {
      create_from_indx(Wdims_c);
    }
      
    get_dims(Wdims);
    assert(Wdims == Wdims_c);

    ecache = new eig_cache();
    assert(ecache);
    assert(!cache->sc);
    cache->sc = ecache;

    ecache->A = new (const double *)[m*m];
    ecache->W = new (double *)[m];
    ecache->A_c = new double[m*m];
    ecache->W_c = new double[m];
    ecache->lwork = 5*m;
    ecache->work = new double[ecache->lwork];

    int jj=0;
    Indx Aindx(2,0),Windx(1,0);
    for (int j=1; j<=m; j++) {
      Windx[0] = j;
      ecache->W[j-1] = location(Windx);
      for (int k=1; k<=m; k++) {
	Aindx[0]=j;
	Aindx[1]=k;
	ecache->A[jj] = A.location(Aindx);
	jj++;
      }
    }
  }

  ecache  = dynamic_cast<eig_cache *> (cache->sc);
  int m = this->dims_p[0].dim;
  const double **pfrom = ecache->A;

  int m2=m*m;
  for (int j=0; j<m2; j++) ecache->A_c[j] = *(ecache->A[j]);
  int info;
  dsyev_("V","U",&m,ecache->A_c,&m,ecache->W_c,
	 ecache->work,&ecache->lwork,&info);

  for (int j=0; j<m; j++) *(ecache->W[j]) = ecache->W_c[j];

  if (!use_cache) delete cache;
  return *this;
#undef A
#undef iA
}
//EOF
_//>

//< print template_subst($seig); //>//

//<=$warn_dont_modify //>
