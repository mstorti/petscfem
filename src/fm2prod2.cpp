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

  // If we concatenate the indices ix0=[indxA,indxB,indxC]
  // we get a vector of size niA+niB+niC. We reorder this
  // in the following way: ix1=[ixfA,ixfB,ixfC,ixcA,ixcB]
  // We have the maps map0to1, and map1to0
  // All of them are in 0 BASE!!
  Indx map0to1(ntot,-1),map1to0(ntot,-1);
  int kk=0;
  for (int j=0; j<niA; j++) {
    int k=ixa[j];
    if (k>0) {
      map1to0[kk] = j;
      map1to0[nfree+kk] = niA+niB+k-1;
      kk++;
    } else map1to0[nfree+niC-1-k] = j;
  }
  map1to0.print("map1to0 after scan A: ");
  for (int j=0; j<niB; j++) {
    int k=ixb[j];
    if (k>0) {
      map1to0[kk] = niA+j;
      map1to0[nfree+kk] = niA+niB+k-1;
      kk++;
    } else map1to0[nfree+niC+nc-1-k] = niA+j;
  }
  map1to0.print("map1to0 after scan B: ");
  // for (int j=0; j<niC; j++) {
  //   int k = Cfdims[j];
  //   map1to0[nfree+k-1] = niA+niB+j;
  // }
  // map1to0.print("map1to0 after scan C: ");
  
  // Compute inverse map
  for (int j=0; j<ntot; j++) 
    map0to1[map1to0[j]] = j;
  map1to0.print("map1to0");
  map0to1.print("map0to1");

  //# Current line ===========       
  exit(0);

#if 0
  // ndimsf.print("dimensions of the free part: ");

  // Dimension C if necessary
  if (!defined) create_from_indx(ndimsf);
  // If the resulting matrix is 2-size then do nothing
  if (comp_storage_size(ndimsf)==0) {
    cache->nelems=0;
    return *this;
  }

  if (ndimsf != Cfdims) {
    Afdims.print("A free dims: ");
    Bfdims.print("B free dims: ");
    ndimsf.print("Combined free dims: ");
    Cfdims.print("Free dims on result matrix: ");
    PETSCFEM_ERROR0("Combined free dims doesn't match free"
                    " dims of  result.\n");
  }

  // Check that the paired contracted indices are
  // equal.
  for (int j=0; j<nc; j++) {
    int k1 = icontr[2*j];
    int k2 = icontr[2*j+1];
    // inA1:= means whether that index
    // belongs to A or to B. Contracted indices
    // must not belong to the same matrix. 
    int nd1,nd2,inA1=0,inA2=0;
    if (k1<=niA) {
      inA1 = 1;
      nd1 = Afdims[k1-1];
    } else {
      nd1 = Bfdims[k1-niA];
    }
    if (k2<=niA) {
      inA2 = 1;
      nd2 = Afdims[k2-1];
    } else {
      nd2 = Bfdims[k2-niA-1];
    }
    ncA *= nd1;
    ncB *= nd1;
    PETSCFEM_ASSERT(nd1==nd2,
                    "Contracted indices must have the same "
                    "dimensions. Offending indices in positions "
                    "%d, and %d.",k1,k2);
    PETSCFEM_ASSERT(inA1!=inA2,
                    "Contracted indices must not be on the "
                    "same matrix. Offending indices in positions "
                    "%d, and %d.",k1,k2);
    ndimsc[j]=nd1;
  }
  printf("nfA %d, ncA %d, nfB %d, ncB %d\n",nfA,ncA,nfB,ncB);

  Indx afindx(nfA,1),acindx(ncA,1);
  while (1) {
    
    if (!inc(afindx,ndimsc)) break;
  }

  // ndimsc.print("dimensions of the contracted part: ");

  Indx findx(nfree,1),cindx(nc,1),tot_indx(ndims,0),
    aindx(niA,0),bindx(niB,0);

  // Loading addresses in cache. For each element in the
  // distination target, we store the complete list of
  // addresses of the lines of elements that contribute to
  // it.
  cache->nlines = mem_size(ndimsf);
  cache->prod_cache.resize(cache->nlines);
  cache->line_cache_start = &*cache->prod_cache.begin();
  cache->line_size = mem_size(ndimsc);
  int line_size = cache->line_size;

    // T 2%
    int jlc=0,inca=1,incb=1;
    LineCache *lc;
    while (1) {
      lc = cache->line_cache_start + jlc++;
      lc->linea.resize(line_size);
      lc->lineb.resize(line_size);
      lc->target = location(findx);

      cindx= Indx(nc,1);
      for (int j=0; j<nfree; j++)
	tot_indx[ifree[j]-1] = findx[j];

      int jj=0;
      lc->linear=0;
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
      lc->starta = &*lc->linea.begin();
      lc->startb = &*lc->lineb.begin();
      lc->inca = inca;
      lc->incb = incb;
      // lc->linear = 0; // force non-linear

      if (!inc(findx,ndimsf)) break;
    }
    // T 97%

    // Hay que contar mejor cuantos elementos hay que traer
    // ctx->op_count.get += cache->nlines*cache->line_size;
    int ntot = cache->nlines*cache->line_size;
    ctx->op_count.put += cache->nlines*cache->line_size;
    ctx->op_count.sum += ntot;
    ctx->op_count.mult += ntot;

    psc = new prod_subcache_t(cache);
    assert(psc);
    assert(!cache->sc);
    cache->sc = psc;
    if (FASTMAT2_USE_DGEMM) {
      psc->ident();
      // if (!psc->not_superlinear_ok()) 
      //   printf("NOT SL!!\n");
      // psc->print();
      
      glob_fm2stats.was_sl_count += psc->superlinear;
      glob_fm2stats.was_not_sl_count += !psc->superlinear;
      
      if (psc->superlinear) {
        glob_fm2stats.was_sl= 1;
        psc->superlinear = glob_fm2stats.use_dgemm;
      }
    }
  }
  fastmat_stats.tcall_not_cached += MPI_Wtime()-start;

  psc = dynamic_cast<prod_subcache_t *> (cache->sc);
  assert(psc);

  if (psc->superlinear) psc->dgemm();
  else {
    // Perform computations using cached addresses
    int 
      nlines = cache->nlines,
      mm=cache->line_size, inca, incb;
    double **pa,**pb,**pa_end,sum,*paa,*pbb,*paa_end;
    LineCache *lc=NULL;

    for (int j=0; j<nlines; j++) {
      lc = cache->line_cache_start+j;
      pa = lc->starta;
      pb = lc->startb;
      inca = lc->inca;
      incb = lc->incb;
      if (lc->linear) {
#if 1
        sum=0.;
        paa = *pa;
        pbb = *pb;
        //         if (inca==1 && incb==1) {
        //           for (int k=0; k<mm; k++) {
        //             sum += (*paa)*(*pbb);
        //             paa++; pbb++;
        //           }
        //         } else 
        if (inca==1) {
          paa_end = paa + mm;
          while (paa<paa_end) {
            sum += (*paa)*(*pbb);
            paa++;
            pbb += incb;
          }
        } else if (incb==1) {
          paa_end = paa + inca*mm;
          while (paa<paa_end) {
            sum += (*paa)*(*pbb);
            paa += inca;
            pbb++;
          }
        } else {
          paa_end = paa + inca*mm;
          while (paa<paa_end) {
            sum += (*paa)*(*pbb);
            paa += inca;
            pbb += incb;
          }
        }
#else
        paa = *pa;
        pbb = *pb;
        paa_end = paa + lc->inca * cache->line_size;
        sum=0.;
        while (paa<paa_end) {
          sum += (*paa)*(*pbb);
          paa += lc->inca;
          pbb += lc->incb;
        }
#endif
      } else {
        pa_end = pa + cache->line_size;
        sum=0.;
        while (pa<pa_end) {
          sum += (**pa++)*(**pb++);
        }
      }
      *(lc->target) = sum;
    }
  }

  if (!ctx->use_cache) delete cache;
  fastmat_stats.tcall += MPI_Wtime()-start;
#endif
  return *this;
}
