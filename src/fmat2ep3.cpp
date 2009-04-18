//<require "prep.pl"//>//
//<=$warn_dont_modify //>

//__INSERT_LICENSE__
//$Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <mkl_cblas.h>

#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fastlib2.h>

int FASTMAT2_USE_DGEMM=1;
int FASTMAT2_USE_DGEMM_VRBS=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class prod_subcache : public FastMatSubCache {
public:
  int superlinear, lda, ldb, ldc, nra, nca, ncb;
  CBLAS_TRANSPOSE transa, transb;
  double *paa0, *pbb0, *pcc0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Helper class, masks the relation betewen index
// and start of `LineCache' or `target' for the function
// `is_superl()' below. 
class superl_helper_t {
public:
  enum  mode_t { a,b,c,none };
  FastMatCache *cache;
  LineCache *line_cache_start;
  superl_helper_t(FastMatCache *cache_a)
    : cache(cache_a), 
      line_cache_start(cache->line_cache_start) { }
  double *address(int j, mode_t mode) {
    if (mode==a) 
      return *line_cache_start[j].starta;
    else if (mode==b) 
      return *line_cache_start[j].startb;
    else if (mode==c) 
      return line_cache_start[j].target;
    else assert(0);
    return NULL;
  }
  int has_rmo(int &m,int &n,int &incrow, int &inccol);
#define MODE_SHORT 0
#define MODE_LONG 1
  int has_rmo2(mode_t mode,
               int &nrow,int &ncol,
               int &incrow, int &inccol,
               int &byrows, int &trvmode);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int superl_helper_t
::has_rmo(int &m,int &n,int &incrow, int &inccol) {

  LineCache *line_cache_start=NULL;
  int N = cache->nlines;
  if (N<1) return 0;
  
  mode_t mode=c;
  double *aa00 = address(0,mode);
  inccol = 1;
  if (N>1) inccol = address(1,mode) - aa00; 
  int dp, j;
  for (j=2; j<N; j++) {
    dp = address(j,mode) - aa00;
    if (dp != j*inccol) break;
  }
  n=j;

  m = N/n;
  incrow = dp;
  if (m==1) incrow = !dp;
  int l = 0;
  for (int j=0; j<m; j++) {
    for (int k=0; k<n; k++) {
      dp = address(l,mode) - aa00; 
      if (dp != incrow*j + inccol*k)
        return 0;
      l++;
    }
  }

  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Determines whether a given set of addresses &a(j), with j
// in[0,N) are in RMO (Rectangular Matrix Order). i.e. there
// exists m,n such that if j = k*n+l and
//
//   &aa(k,l) = &aa(0,0) + inccol*k + rowcol*l
int superl_helper_t
::has_rmo2(mode_t mode,int &nrow,int &ncol,
           int &incrow, int &inccol,
           int &byrows, int &trvmode) {

  int inc1=0, inc2=0, size1, size2;
  LineCache *line_cache_start=NULL;
  int N = cache->nlines;

  double *aa00 = address(0,mode);
  inc1 = 1;
  if (N>1) inc1 = address(1,mode) - aa00; 
  inccol = (mode==a? cache->line_cache_start->inca :
            mode==b? cache->line_cache_start->incb : -1);
  assert(inccol!=-1);

  ncol = cache->line_size;

  int dp, j;
  for (j=2; j<N; j++) {
    dp = address(j,mode) - aa00;
    if (dp != j*inc1) break;
  }

  size2=j;
  if (N%size2 != 0) return 0;
  size1 = N/size2;
  inc2 = dp;
  if (size1==1) inc2 = !dp;
  int l = 0;
  for (int j=0; j<size1; j++) {
    for (int k=0; k<size2; k++) {
      dp = address(l,mode) - aa00; 
      if (dp != inc2*j + inc1*k)
        return 0;
      l++;
    }
  }
  
  // It is a matrix
  assert(inc1==0 || inc2==0);
  
  byrows = 1;
  if (inc1==0) {
    trvmode = MODE_SHORT;
    incrow = inc2;
    nrow = size1;
  } else {
    trvmode = MODE_LONG;
    incrow = inc1;
    nrow = size2;
  }

  if (inccol!=1) {
    if (incrow!=1) return 0;
    else {
      // swap rows with cols
      swap(inccol,incrow);
      swap(ncol,nrow);
      byrows = !byrows;
    }
  }

#if 0  
  printf("nrow %d, ncol %d, incrow %d, inccol %d, trvmode %s, "
         "byrows %d\n",nrow,ncol,incrow,inccol,
         (trvmode==MODE_SHORT? "short" : "long"),byrows);
#endif

  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
//<$prod = <<'//EOF';
FastMat2 & FastMat2::prod(const FastMat2 & A,const FastMat2 & B,
                          const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("__NAME__",this,&A,&B);
    ctx->check(m);
    Indx indx;
    READ_ARG_LIST(arg,indx,INT_ARG_LIST_DEFAULT_VAL,EXIT2)
    ctx->check(indx);
  }
#endif
  FastMatCache *cache = ctx->step();

  prod_subcache *psc=NULL;
  if (!ctx->was_cached  ) {
    Indx ia,ib,ii;

    Indx Afdims,Bfdims,fdims;
    A.get_dims(Afdims);
    B.get_dims(Bfdims);

    // maxc:= maximum contracted index
    ii.push_back(m);
    int niA = Afdims.size();
    int niB = Bfdims.size();
    int ndims = niA+niB;

#ifdef USE_VAR_ARGS
    va_start(ap,m);
    read_int_list(niA+niB-1,ap,&ii);
    va_end(ap);
#else
    READ_INT_ARG_LIST(ii);
    assert(ii.size() == niA+niB);
#endif

    int nfree=0,nc=0;
    for (int j=0; j<ii.size(); j++) {
      int k = ii[j];
      if (k>0 && k>nfree) nfree = k;
      if (k<0 && -k>nc) nc = -k;
    }

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
  
    Indx ndimsf(nfree,0),ndimsc(nc,0);
    for (int j=0; j<nfree; j++) {
      int k = ifree[j];
      if (k<=niA) {
	ndimsf[j] = Afdims[k-1];
      } else {
	ndimsf[j] = Bfdims[k-niA-1];
      }
    }

    // ndimsf.print("dimensions of the free part: ");

  // Dimension C if necessary
    if (!defined) create_from_indx(ndimsf);
    if (comp_storage_size(ndimsf)==0) {
      cache->nelems=0;
      return *this;
    }

    get_dims(fdims);
    if (ndimsf != fdims) {
      Afdims.print("A free dims: ");
      Bfdims.print("B free dims: ");
      ndimsf.print("Combined free dims: ");
      fdims.print("Free dims on result matrix: ");
      printf("Combined free dims doesn't match free"
	     " dims of  result.\n");
      abort();
    }
    for (int j=0; j<nc; j++) {
      int k1 = icontr[2*j];
      int k2 = icontr[2*j+1];
      int nd1,nd2;
      if (k1<=niA) {
	nd1 = Afdims[k1-1];
      } else {
	nd1 = Bfdims[k1-niA];
      }
      if (k2<=niA) {
	nd2 = Afdims[k2-1];
      } else {
	nd2 = Bfdims[k2-niA-1];
      }
      assert(nd1==nd2);
      ndimsc[j]=nd1;
    }

    // ndimsc.print("dimensions of the contracted part: ");

    Indx findx(nfree,1),cindx(nc,1),tot_indx(ndims,0),
      aindx(niA,0),bindx(niB,0);

  // Loading addresses in cache
  // For each element in the distination target, we store the complete
  // list of addresses of the lines of elements that contribute to
  // it. 
    cache->nlines = mem_size(ndimsf);
    cache->prod_cache.resize(cache->nlines);
    cache->line_cache_start = &*cache->prod_cache.begin();
    cache->line_size = mem_size(ndimsc);
    int line_size = cache->line_size;

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
    // Hay que contar mejor cuantos elementos hay que traer
    // ctx->op_count.get += cache->nlines*cache->line_size;
    int ntot = cache->nlines*cache->line_size;
    ctx->op_count.put += cache->nlines*cache->line_size;
    ctx->op_count.sum += ntot;
    ctx->op_count.mult += ntot;

    int superlinear = 0, lda=-1, ldb=-1, ldc, 
      nra=-1, nca=-1, nrb=-1;
    double *paa0=NULL, *pbb0=NULL, *pcc0=NULL;
    {
      int sl, 
        nrowa,ncola,incrowa,inccola,byrowsa,trvmodea,
        nrowb,ncolb,incrowb,inccolb,byrowsb,trvmodeb,
        nrowc,ncolc,incrowc,inccolc,byrowsc,trvmodec;
      int nrowopa,ncolopa,nrowopb,ncolopb,transa,transb;
      LineCache *lc0 = cache->line_cache_start;
      double *paa0,*pbb0,*pcc0;
        
      superl_helper_t obj(cache);

      sl = obj.has_rmo2(superl_helper_t::a,
                        nrowa,ncola,incrowa,inccola,byrowsa,trvmodea);
      if (!sl) goto NOT_SL;
      sl = obj.has_rmo2(superl_helper_t::b,
                        nrowb,ncolb,incrowb,inccolb,byrowsb,trvmodeb);
      if (!sl) goto NOT_SL;
      
      // opa: refers to op(A) in BLAS documentation
      nrowopa=nrowa; ncolopa=ncola;
      nrowopb=nrowb; ncolopb=ncolb;
      transa=0; transb=0;
      
      if (!byrowsa) {
        swap(nrowopa,ncolopa);
        transa=1;
      }
      if (byrowsb) {
        swap(nrowopb,ncolopb);
        transb=1;
      }
      assert(ncolopa==nrowopb);
#if 0
      printf("args to dgemm:\n");
      printf("   A: %d x %d, lda %d, trans %d\n",
             nrowopa,ncolopa,incrowa,transa);
      printf("   B: %d x %d, lda %d, trans %d\n",
             nrowopb,ncolopb,incrowb,transb);
#endif
      sl = obj.has_rmo(nrowc,ncolc,incrowc,inccolc);
      // If C is traversed by columns then we should transpose all the
      // product: C' = B'*A', e.g. transpose A, B and exchange
      if (!sl) goto NOT_SL;
      assert(nrowc*ncolc==nrowopa*ncolopb);
      if (nrowc!=nrowopa && nrowc==1) {
        nrowc=nrowopa;
        ncolc=ncolopb;
        incrowc=inccolc*ncolc;
      }
      if (ncolc==1) inccolc=1;

      superlinear = FASTMAT2_USE_DGEMM;

      lc0 = cache->line_cache_start;
      paa0 = *lc0->starta;
      pbb0 = *lc0->startb;
      pcc0 = lc0->target;

    NOT_SL: ;
    
      if (FASTMAT2_USE_DGEMM && FASTMAT2_USE_DGEMM_VRBS) 
        printf("use dgemm (superlinear) %d\n",superlinear);

      psc = new prod_subcache;
      assert(psc);
      assert(!cache->sc);
      cache->sc = psc;
      psc->superlinear = superlinear;
      if (superlinear) {
        psc->nra = nrowopa;
        psc->nca = ncolopa;
        psc->lda = incrowa;
        psc->ncb = ncolopb;
        psc->ldb = incrowb;
        psc->ldc = incrowc;
        psc->paa0 = paa0;
        psc->pbb0 = pbb0;
        psc->pcc0 = pcc0;
        psc->transa = (transa? CblasTrans : CblasNoTrans);
        psc->transb = (transb? CblasTrans : CblasNoTrans);
      }
    }
  }

  psc = dynamic_cast<prod_subcache *> (cache->sc);
  assert(psc);

  if (psc->superlinear) {
    int p = cache->line_size;
    cblas_dgemm(CblasRowMajor,psc->transa,psc->transb,
                psc->nra,psc->ncb,psc->nca,1.0,
                psc->paa0,psc->lda,psc->pbb0,psc->ldb,0.0,
                psc->pcc0,psc->ldc);
  } else {
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
  return *this;
}
//EOF
_//>

//< print template_subst($prod); //>//

//<=$warn_dont_modify //>
