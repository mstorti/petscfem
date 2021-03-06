// DON'T EDIT MANUALLY THIS FILE !!!!!!
// This files automatically generated by ePerl from 
// the corresponding file (omitting the _eperl) name.


//__INSERT_LICENSE__
//$Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$
#include <math.h>
#include <stdio.h>

#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fastlib2.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class ns_eig_cache : public FastMatSubCache {
public:
  const double **A;
  double **W, **VR, **VL;
  double *A_c,*work,*W_c, *VR_c, *VL_c;
  int lwork, m;
  int crev, clev;
  ns_eig_cache() { 
    W=VR=VL=NULL; A=NULL; 
    A_c=W_c=VR_c=VL_c=NULL; work=NULL; 
    crev=0; clev=0;
  }
  ~ns_eig_cache() { 
    delete[] A;
    delete[] W;
    delete[] VR;
    delete[] VL;
    delete[] A_c;
    delete[] W_c;
    delete[] VR_c;
    delete[] VL_c;
    delete[] work;
  }
};

extern "C"
void dgeev_(const char *jobvl,const char *jobvr,int *n,double *a,int *lda, 
	    double *wr,double *wi,double *vl,int *ldvl,double *vr,int *ldvr,
	    double *work,int *lwork,int *info);



//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 &
FastMat2::eig(const FastMat2 & A, 
	      FastMat2 *VR, FastMat2 *VL) { 

  CTX2_CHECK("eig",this,&A,VR);

  ns_eig_cache *sc;
  int symmetric = 0;
  if (!ctx->was_cached) {
    Indx Adims;
    A.get_dims(Adims);
    int ndims = Adims.size();
    assert (ndims==2);
    assert (Adims[0] == Adims[1]);
    int m= Adims[0];

    Indx Wdims, Wdims_c, Vdims;
    Wdims_c.push_back(2);
    Wdims_c.push_back(m);
    if (!defined) create_from_indx(Wdims_c);
      
    get_dims(Wdims);
    assert(Wdims == Wdims_c);
  
    sc = new ns_eig_cache();
    assert(sc);
    assert(!cache->sc);
    cache->sc = sc;

    if (VL) sc->clev=1;
    if (sc->clev) {
      assert(VL);
      if (!VL->defined) VL->create_from_indx(Adims);
      VL->get_dims(Vdims);
      assert(Adims == Vdims);
    }

    if (VR) sc->crev=1;
    if (sc->crev) {
      assert(VR);
      if (!VR->defined) VR->create_from_indx(Adims);
      VR->get_dims(Vdims);
      assert(Vdims == Adims);
    }

    sc->m = m;
    sc->A = new const double *[m*m];
    sc->W = new double *[2*m];
    if (sc->clev) {
      sc->VL = new double *[m*m];
      sc->VL_c = new double[m*m];
    }
    if (sc->crev) {
      sc->VR = new double *[m*m];
      sc->VR_c = new double[m*m];
    }
    sc->A_c = new double[m*m];
    sc->W_c = new double[2*m];
    sc->lwork = 5*m;
    sc->work = new double[sc->lwork];

    int jj=0;
    Indx Aindx(2,0),Windx(2,0);
    for (int j=1; j<=m; j++) {
      Windx[1] = j;
      Windx[0] = 1;
      sc->W[j-1] = location(Windx);
      Windx[0] = 2;
      sc->W[m+j-1] = location(Windx);
      for (int k=1; k<=m; k++) {
	Aindx[0]=k;
	Aindx[1]=j;
	sc->A[jj] = A.location(Aindx);
	if (sc->clev) sc->VL[jj] = VL->location(Aindx);
	if (sc->crev) sc->VR[jj] = VR->location(Aindx);
	jj++;
      }
    }
  }

  sc  = dynamic_cast<ns_eig_cache *> (cache->sc);
  int &m = sc->m;
  const double **pfrom = sc->A;

  int m2=m*m;
  for (int j=0; j<m2; j++) sc->A_c[j] = *(sc->A[j]);
  int info;
  
  dgeev_((sc->clev ? "V" : "N"), (sc->crev ? "V" : "N"),
	 &m,sc->A_c, &m, 
	 sc->W_c,sc->W_c+m,sc->VL_c,&m,sc->VR_c,&m,
	 sc->work,&sc->lwork,&info);
  
  for (int j=0; j<2*m; j++) *(sc->W[j]) = sc->W_c[j];
  if (sc->clev) for (int j=0; j<m2; j++) *(sc->VL[j]) = sc->VL_c[j];
  if (sc->crev) for (int j=0; j<m2; j++) *(sc->VR[j]) = sc->VR_c[j];

  if (!ctx->use_cache) delete cache;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::eig(const FastMat2 & A, FastMat2 &VR) {
  eig(A,&VR);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::eig(const FastMat2 & A, FastMat2 &VR, FastMat2 &VL) {
  eig(A,&VR,&VL);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class eig_cache : public FastMatSubCache {
public:
  const double **A;
  double **W, **V;
  double *A_c,*work,*W_c;
  int lwork;
  eig_cache() { W=V=NULL; A=NULL; A_c=W_c=NULL; work=NULL; }
  ~eig_cache() { 
    delete[] A;
    delete[] W;
    delete[] V;
    delete[] A_c;
    delete[] W_c;
    delete[] work;
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & 
FastMat2::seig(const FastMat2 & A) {
  seig(A,(FastMat2 & )A,0);
  return *this; 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
extern "C" void dsyev_(const char*jobz,const char *uplo,int *n,double *a,int *lda,
		       double *w,double *work, int *lwork, int *info);



//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 &
FastMat2::seig(const FastMat2 & A, FastMat2 &V,
               int compute_eigen_vectors) { 

  CTX2_CHECK("seig",this,&A,&V);

  eig_cache *ecache;
  int &cev = compute_eigen_vectors;
  if (!ctx->was_cached) {
    Indx Adims;
    A.get_dims(Adims);
    int ndims = Adims.size();
    assert (ndims==2);
    assert (Adims[0] == Adims[1]);
    int m= Adims[0];

    Indx Wdims,Wdims_c, Vdims;
    Wdims_c.push_back(m);
    if (!defined) {
      create_from_indx(Wdims_c);
    }
      
    get_dims(Wdims);
    assert(Wdims == Wdims_c);
    
    if (cev) {
      if (!V.defined) V.create_from_indx(Adims);
      V.get_dims(Vdims);
      assert(Adims == Vdims);
    }

    ecache = new eig_cache();
    assert(ecache);
    assert(!cache->sc);
    cache->sc = ecache;

    ecache->A = new const double *[m*m];
    ecache->W = new double *[m];
    if (cev) ecache->V = new double *[m*m];
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
	Aindx[1]=j;
	Aindx[0]=k;
	ecache->A[jj] = A.location(Aindx);
	if (cev) ecache->V[jj] = V.location(Aindx);
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
  dsyev_((cev ? "V" : "N"),
	 "U",&m,ecache->A_c,&m,ecache->W_c,
	 ecache->work,&ecache->lwork,&info);

  for (int j=0; j<m; j++) *(ecache->W[j]) = ecache->W_c[j];
  if (cev) for (int j=0; j<m2; j++) *(ecache->V[j]) = ecache->A_c[j];

  if (!ctx->use_cache) delete cache;
  return *this;
}

class inv_cache : public FastMatSubCache {
public:
  double *A,*b,*inv_A;
  int *ipvt;
  inv_cache() : A(NULL), b(NULL), inv_A(NULL), ipvt(NULL) {}
  ~inv_cache() {
    delete[] A;
    delete[] inv_A;
    delete[] ipvt;
    delete[] b;
  }
};

extern "C" void dgefa_(double *a,int *lda,int *n,int *ipvt,int *info);
extern "C" void dgesl_(double *a,int *lda,int *n,int *ipvt,double *b,int *job);


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::inv(const FastMat2 & A) {

  CTX2_CHECK("inv",this,&A);

  if (!ctx->was_cached) {
    Indx Adims;
    A.get_dims(Adims);
    int ndims = Adims.size();
    assert (ndims==2);
    assert (Adims[0] == Adims[1]);
    int m= Adims[0];

    if (!defined)
      create_from_indx(Adims);
      
    Indx dims_;
    get_dims(dims_);
    assert(dims_ == Adims);

    if (m<=3) {
      cache->nelems=m*m;
      cache->from_elems.resize(cache->nelems);
      cache->pfrom = &*cache->from_elems.begin();
      cache->to_elems.resize(cache->nelems);
      cache->pto = &*cache->to_elems.begin();
      int jj=0;
      Indx indx(2,0);
      for (int j=1; j<=m; j++) {
	for (int k=1; k<=m; k++) {
	  indx[0]=j;
	  indx[1]=k;
	  cache->pto[jj] = location(indx);
	  cache->pfrom[jj] = A.location(indx);
	  jj++;
	}
      }
    } else {
#define USE_LAPACK
#ifndef USE_LAPACK
      cache->A = new Matrix(m,m);
      cache->B = new Matrix(m,m);
#else
      inv_cache *isc = new inv_cache;
      assert(isc);
      isc->A = new double[m*m];
      isc->inv_A = new double[m*m];
      isc->b = new double[m];
      isc->ipvt = new int[2*m];
      cache->sc = isc;
#endif
    }

    // This I don't know if it's correct.
    ctx->op_count.get += m*m;
    ctx->op_count.put += m*m;
    ctx->op_count.mult += m*m*m;
    ctx->op_count.sum += m*m*m;

  }

  int m = this->dims_p[0].dim;
  double det;
  double **pfrom = cache->pfrom;
  double **pto = cache->pto;
#define A(i,j) (*pfrom[(i-1)*M+(j-1)])
#define iA(i,j) (*pto[(i-1)*M+(j-1)])
  if (m==1) {
    **pto = 1. / **pfrom;
  } else if (m==2) {
#define M 2
    det = A(1,1)*A(2,2)-A(1,2)*A(2,1);
    iA(1,1) = A(2,2)/det;
    iA(1,2) = -A(1,2)/det;
    iA(2,1) = -A(2,1)/det;
    iA(2,2) = A(1,1)/det;
#undef M
  } else if (m==3) {
#define M 3
    double c11,c21,c31;

    c11 = A(2,2)*A(3,3)-A(2,3)*A(3,2);
    c21 = A(2,3)*A(3,1)-A(2,1)*A(3,3);
    c31 = A(2,1)*A(3,2)-A(2,2)*A(3,1);

    det= A(1,1) * c11 + A(1,2) * c21 + A(1,3) * c31;

    iA(1,1) = c11/det;
    iA(1,2) = (A(3,2)*A(1,3)-A(3,3)*A(1,2))/det;
    iA(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))/det;
	   
    iA(2,1) = c21/det;
    iA(2,2) = (A(3,3)*A(1,1)-A(3,1)*A(1,3))/det;
    iA(2,3) = (A(1,3)*A(2,1)-A(1,1)*A(2,3))/det;
	   
    iA(3,1) = c31/det;
    iA(3,2) = (A(3,1)*A(1,2)-A(3,2)*A(1,1))/det;
    iA(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/det;
#undef M
  } else if (m>3) {
#ifndef USE_LAPACK
    A.export_vals(*cache->A);
    *cache->B = cache->A->i();
    set(cache->B->Store());
#else
    inv_cache *isc = dynamic_cast<inv_cache *>(cache->sc);
    A.export_vals(isc->A);
    int info;
    dgefa_(isc->A,&m,&m,isc->ipvt,&info);
    assert(!info);
    int job = 1;
    for (int k=0; k<m; k++) {
      for (int l=0; l<m; l++) isc->b[l] = 0.;
      isc->b[k] = 1.;
      dgesl_(isc->A,&m,&m,isc->ipvt,isc->b,&job);
      for (int l=0; l<m; l++) isc->inv_A[m*l+k] = isc->b[l];
    }
    set(isc->inv_A);
#endif
  }

  if (!ctx->use_cache) delete cache;
  return *this;
#undef A
#undef iA
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double FastMat2::trace() {

  CTX2_CHECK("trace",this);

  if (!ctx->was_cached) {
    assert(defined);
    Indx Adims;
    get_dims(Adims);
    int ndims = Adims.size();
    assert (ndims==2);
    assert (Adims[0] == Adims[1]);
    int n= Adims[0];

    Indx indx;
    cache->to_elems.resize(n);
    for (int j=1; j<=n; j++) {
      indx = Indx(2,j);
      cache->to_elems[j-1] = location(indx);
    } 
    cache->nelems = n;
    cache->pto = &*cache->to_elems.begin();
    ctx->op_count.put += n;
    ;
  
  }

  double **pto = cache->pto;
  double **pto_end = pto + cache->nelems;

  double tr = 0.0;
  while (pto < pto_end) {
    tr += **pto++;
  }

  if (!ctx->use_cache) delete cache;
  return tr;
}

// DON'T EDIT MANUALLY THIS FILE !!!!!!
// This files automatically generated by ePerl from 
// the corresponding file (omitting the _eperl) name.

