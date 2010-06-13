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

#if 0
  // Count how many free and contracted
  // indices are there
  int nfree=0,nc=0;
  for (int j=0; j<ii.size(); j++) {
    int k = ii[j];
    if (k>0 && k>nfree) nfree = k;
    if (k<0 && -k>nc) nc = -k;
  }

  // Nbr of free indices in A and B
  int nfreeA=0,nfreeB=0;
  for (int j=0; j<niA; j++) 
    nfreeA += ii[j]>0;
  for (int j=niA; j<niA; j++) 
    nfreeA += ii[j]>0;

  // Get the free and contracted indices
  Indx ifree(nfree,0),icontr(2*nc,0);
  for (int j=0; j<ii.size(); j++) {
    int k = ii[j];
    if (k>0) {
      ifree[k-1] = j+1;
    } else {
      k = -k;
      if (icontr[2*(k-1)]==0) {
        icontr[2*(k-1)]=j+1;
      } else {
        icontr[2*k-1]=j+1;
      }
    }
  }
  // ifree.print("ifree: ");
  // icontr.print("icontr: ");
    
  // ndimsf:= dimensions of the free indices
  // ndimsc:= dimensions of the contracted indices
  Indx ndimsf(nfree,0),ndimsc(nc,0);
  int nfA=1,ncA=1,nfB=1,ncB=1;
  for (int j=0; j<nfree; j++) {
    int k = ifree[j];
    if (k<=niA) {
      ndimsf[j] = Afdims[k-1];
      nfA *= ndimsf[j];
    } else {
      ndimsf[j] = Bfdims[k-niA-1];
      nfB *= ndimsf[j];
    }
  }
#endif

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
  Indx mapA(niA,-1),mapB(niB,-1),mapC(niC,-1);
  int kf=0;
  for (int j=0; j<niA; j++) {
    int k=ixa[j];
    if (k>0) {
      mapA[kf] = j;
      mapC[kf] = k-1;
      kf++;
    } else mapA[nfA-1-k] = j;
  }
  // printf("After scanning A\n");
  // mapA.print("mapA: ");
  // mapB.print("mapB: ");
  // mapC.print("mapC: ");
  kf=0;
  for (int j=0; j<niB; j++) {
    int k=ixb[j];
    if (k>0) {
      mapB[nc+kf] = j;
      mapC[nfA+kf] = k-1;
      kf++;
    } else mapB[-1-k] = j;
  }
  // printf("After scanning B\n");
  mapA.print("mapA: ");
  mapB.print("mapB: ");
  mapC.print("mapC: ");

#if 0
  while(1) {
    for (int j=0; j<nc; j++) {
      int k1=icontr[2*j];
      int k2=icontr[2*j+1];
      tot_indx[k1-1] = cindx[j];
      tot_indx[k2-1] = cindx[j];
    }

    // Extract the A and B parts of the indices
    copy(&tot_indx[0],&tot_indx[niA],&aindx[0]);
    copy(&tot_indx[niA],&tot_indx[ndims],&bindx[0]);

    lc->linea[jj] = A.location(aindx);
    lc->lineb[jj] = B.location(bindx);

    if (jj==1) {
      inca = lc->linea[1] - lc->linea[0];
      incb = lc->lineb[1] - lc->lineb[0];
      lc->linear=1;
    } else if (lc->linear && jj>1) {

      if (lc->linea[jj] - lc->linea[jj-1] != inca) lc->linear=0;
      if (lc->lineb[jj] - lc->lineb[jj-1] != incb) lc->linear=0;
    }

    jj++;
    if (!inc(cindx,ndimsc)) break;
  }
#endif  
  
  //# Current line ===========       
  exit(0);
  return *this;
}
