//<require "prep.pl"; require "fmat2.pl"; //>//
//<=$warn_dont_modify //>

//__INSERT_LICENSE__
//$Id: fmat2ep2.cpp,v 1.5 2002/11/29 12:56:29 mstorti Exp $
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
  ns_eig_cache() { 
    W=VR=VL=NULL; A=NULL; 
    A_c=W_c=VR_c=VL_c=NULL; work=NULL; 
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

//<$eig=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 &
FastMat2::eig(const FastMat2 & A, FastMat2 *VR=NULL, FastMat2 *VL=NULL, 
	      int crev=0, int clev=0) { 

  __CACHE_OPERATIONS__;

  ns_eig_cache *sc;
  int symmetric = 0;
  if (!was_cached) {
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
  
    if (VL) clev=1;
    if (clev) {
      assert(VL);
      if (!VL->defined) VL->create_from_indx(Adims);
      VL->get_dims(Vdims);
      assert(Adims == Vdims);
    }

    if (VR) crev=1;
    if (crev) {
      assert(VR);
      if (!VR->defined) VR->create_from_indx(Adims);
      VR->get_dims(Vdims);
      assert(Vdims == Adims);
    }

    sc = new ns_eig_cache();
    assert(sc);
    assert(!cache->sc);
    cache->sc = sc;

    sc->m = m;
    sc->A = new (const double *)[m*m];
    sc->W = new (double *)[2*m];
    if (clev) {
      sc->VL = new (double *)[m*m];
      sc->VL_c = new double[m*m];
    }
    if (crev) {
      sc->VR = new (double *)[m*m];
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
	if (clev) sc->VL[jj] = VL->location(Aindx);
	if (crev) sc->VR[jj] = VR->location(Aindx);
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
  
  dgeev_((clev ? "V" : "N"), (crev ? "V" : "N"),
	 &m,sc->A_c, &m, 
	 sc->W_c,sc->W_c+m,sc->VL_c,&m,sc->VR_c,&m,
	 sc->work,&sc->lwork,&info);
  
  for (int j=0; j<2*m; j++) *(sc->W[j]) = sc->W_c[j];
  if (clev) for (int j=0; j<m2; j++) *(sc->VL[j]) = sc->VL_c[j];
  if (crev) for (int j=0; j<m2; j++) *(sc->VR[j]) = sc->VR_c[j];

  if (!use_cache) delete cache;
  return *this;
}
//EOF
_//>

//< print template_subst($eig); //>//

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

//<$seig=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 &
FastMat2::seig(const FastMat2 & A, FastMat2 &V,int compute_eigen_vectors=1) { 

  __CACHE_OPERATIONS__;

  eig_cache *ecache;
  int &cev = compute_eigen_vectors;
  if (!was_cached) {
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

    ecache->A = new (const double *)[m*m];
    ecache->W = new (double *)[m];
    if (cev) ecache->V = new (double *)[m*m];
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

  if (!use_cache) delete cache;
  return *this;
}
//EOF
_//>

//< print template_subst($seig); //>//

//<=$warn_dont_modify //>
