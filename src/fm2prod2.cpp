//__INSERT_LICENSE__
#include <cmath>
#include <ctime>
#include <cstdio>
#include <unistd.h>

#include <algorithm>
#ifdef USE_MKL
#include <mkl_cblas.h>
#endif

#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fm2stats.h>
#include <src/fastlib2.h>
#include <src/fm2prod.h>

int FASTMAT2_USE_PROD2=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class prod2_subcache_t : public FastMatSubCache {
public:
  vector<double *> ap,bp,cp;  
  vector<double> a,b,c;
  int nA,nB,nC, nrowa,ncola,ncolb;
  void make_prod();
  prod2_subcache_t() {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::get_addresses(Indx permA,Indx Afdims,
                             vector<double *> &ap) const {
  int n = permA.size();
  Indx iA(n,1), aindx(n,1),dA(n,-1);
  for (int j=0; j<n; j++) dA[j] = Afdims[permA[j]];
  int j=0;
  while (1) {
    for (int l=0; l<n; l++) {
      int k = permA[l];
      aindx[k] = iA[l];
    }
    ap[j++] = location(aindx);
    if (!inc(iA,dA)) break;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// This is the function that makes the product of two matrices.
// The others for 3,4,etc... are wrappers to this one. 
// Computes C = A*B, where C is *this
FastMat2 & 
FastMat2::prod2(const FastMat2 &A,const FastMat2 &B,
                vector<int> &ixa, 
                vector<int> &ixb) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod2",this,&A,&B);
    Indx Ixa,Ixb;
    for (unsigned int j=0; j<ixa.size(); j++)
      Ixa.push_back(ixa[j]);
    for (unsigned int j=0; j<ixb.size(); j++)
      Ixb.push_back(ixb[j]);
    ctx->check(Ixa);
    ctx->check(Ixb);
  }
#endif
  FastMatCache *cache = ctx->step();

  prod2_subcache_t *psc = NULL;
  if (!ctx->was_cached  ) {

    psc = new prod2_subcache_t;
    assert(psc);
    assert(!cache->sc);
    cache->sc = psc;

    // get the free indices of A and B
    Indx ii,Afdims,Bfdims,Cfdims;
    A.get_dims(Afdims);
    B.get_dims(Bfdims);
  
    int niA = Afdims.size();
    int niB = Bfdims.size();
    int ndims = niA+niB;

    int ixas = int(ixa.size());
    int ixbs = int(ixb.size());
    PETSCFEM_ASSERT(ixas == niA,
                    "Nbr of indices passed in the contraction for A"
                    " must be equal to the actual indices of A. "
                    "nbr of A indices %d. Passed %d",ixas,niA);  
    PETSCFEM_ASSERT(ixbs == niB,
                    "Nbr of indices passed in the contraction for B"
                    " must be equal to the actual indices of B. "
                    "nbr of B indices %d. Passed %d",ixbs,niB);  
    
    // This implementation is based on an older one
    // that used the vector of indices for contraction.
    // Then we simply put the `ixa' and `ixb' indices in `ii'
    // and continue. 
    for (unsigned int j=0; j<ixa.size(); j++) 
      ii.push_back(ixa[j]);
    for (unsigned int j=0; j<ixb.size(); j++) 
      ii.push_back(ixb[j]);

    // Get free dims of output matrix
    get_dims(Cfdims);
    int niC = Cfdims.size();

    // nfree:= nbr of free indices
    // nctr:= nbr of contracted indices
    // Id we put the indices in A,B, and C
    // we should have niA+niB+niC = 2*(nfree+nc)
    // Compute permutation of indices
    // indx_perm[old_pos] -> new_pos
    // So that in the new position we have
    int nfA=0,nfB=0,nfree=0,nc=0;
    for (int j=0; j<niA; j++) {
      if (ixa[j]>0) nfA++;
      else nc++;
    }
    for (int j=0; j<niB; j++) {
      if (ixb[j]>0) nfB++;
      else nc++;
    }
    PETSCFEM_ASSERT(nc%2==0,
                    "Nbr of contracted indices must be even."
                    "nc  %d",nc);  

    nc /= 2;
    nfree = nfA+nfB;
    int ntot = niA+niB+niC;
    PETSCFEM_ASSERT(ntot==2*(nfree+nc),
                    "Bad balance of indices. "
                    "niA %d, niB %d, niC %d, nfree %d, nctr %d",
                    niA,niB,niC,nfree,nc);  
  
    // We reorder the indices of A in (fA,cA), B in (cB,fB). 
    // `f' denotes free and `c' contracted.
    // Indices in fA and fB keep the order as was given.
    // cA and cB are reordered according to the negative
    // index that was given. The indices in C are also
    // reordered so as to keep the same order as appear
    // in fA,fB concatenated.
    // Example: if we make C(lim)=A(kij)B(lkmj) 
    // then: we reorder this as
    // C(ilm)=A(ijk)B(kjlm) 
    // So that we have the following permutations
    // mapC=[1,0,2] mapA=[2,0,1]
    int nrowa=1,ncola=1,nrowb=1,ncolb=1;
    Indx mapA(niA,-1),mapB(niB,-1),mapC(niC,-1);
    int kf=0;
    for (int j=0; j<niA; j++) {
      int k=ixa[j];
      if (k>0) {
        mapA[kf] = j;
        mapC[kf] = k-1;
        kf++;
        nrowa *= Afdims[j];
      } else {
        mapA[nfA-1-k] = j;
        ncola *= Afdims[j];
      }
    }

    kf=0;
    for (int j=0; j<niB; j++) {
      int k=ixb[j];
      if (k>0) {
        mapB[nc+kf] = j;
        mapC[nfA+kf] = k-1;
        kf++;
        ncolb *= Bfdims[j];
      } else {
        mapB[-1-k] = j;
        nrowb *= Bfdims[j];
      }
    }
    assert(ncola==nrowb);

    int 
      nA = A.size(),
      nB = B.size(),
      nC = size();

    vector<double *> 
      &ap = psc->ap,
      &bp = psc->bp,
      &cp = psc->cp;

    vector<double> 
      &a = psc->a,
      &b = psc->b,
      &c = psc->c;

    ap.resize(nA);
    bp.resize(nB);
    cp.resize(nC);

    A.get_addresses(mapA,Afdims,ap);
    B.get_addresses(mapB,Bfdims,bp);
    get_addresses(mapC,Cfdims,cp);

    a.resize(nA,0.0);
    b.resize(nB,0.0);
    c.resize(nC,0.0);

    psc->nA = nA;
    psc->nB = nB;
    psc->nC = nC;

    psc->nrowa = nrowa;
    psc->ncola = ncola;
    psc->ncolb = ncolb;
  }

  psc = dynamic_cast<prod2_subcache_t *>(cache->sc);
  assert(psc);

  psc->make_prod();
  
  if (!ctx->use_cache) delete cache;
  return *this;
}

void prod2_subcache_t::make_prod() {
  for (int j=0; j<nA; j++) a[j] = *ap[j];
  for (int j=0; j<nB; j++) b[j] = *bp[j];
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
              nrowa,ncolb,ncola,1.0,&a[0],ncola,&b[0],ncolb,0.0,
              &c[0],ncolb);
  for (int j=0; j<nC; j++) *cp[j] = c[j];
}
