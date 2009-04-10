//<require "prep.pl"//>//
//<=$warn_dont_modify //>

//__INSERT_LICENSE__
//$Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$
#include <math.h>
#include <stdio.h>
#include <mkl_cblas.h>

#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fastlib2.h>

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
  enum  mode_t { a,b,c };
  LineCache *line_cache_start;
  superl_helper_t(LineCache *line_cache_start_a)
    : line_cache_start(line_cache_start_a) { }
  double *address(int j, mode_t mode) {
    if (mode==a) 
      return *line_cache_start[j].starta;
    else if (mode==b) 
      return *line_cache_start[j].startb;
    else if (mode==c) 
      return line_cache_start[j].target;
    else assert(0);
  }
  int has_rmo(int N,mode_t mode,int &nrow,int &ncol, 
              int &inccol, int &incrow);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Determines whether a given set of addresses &a(j), with j
// in[0,N) are in RMO (Rectangular Matrix Order). i.e. there
// exists m,n such that if j = k*n+l and
//
//   &aa(k,l) = &aa(0,0) + inccol*k + rowcol*l
int superl_helper_t
::has_rmo(int N,mode_t mode,int &m,int &n, 
          int &incrow, int &inccol) {
  m=1; n=2; inccol=1; incrow=1;
  if (N<2) return 1;

  double *aa00 = address(0,mode);
  inccol = address(1,mode) - aa00; 
  int dp, j;
  m=1; incrow = N;
  for (j=2; j<N; j++) {
    dp = address(j,mode) - aa00;
    if (dp != j*inccol) break;
  }
  n=j;
  if (j<N) {
    if (N%n!=0) return 0;
    m = N/n;
    incrow = dp;
    int l = 0;
    for (int j=0; j<m; j++) {
      for (int k=0; k<n; k++) {
        dp = address(l,mode) - aa00; 
        if (dp != incrow*j + inccol*k)
          return 0;
        l++;
      }
    }
  }
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
      superl_helper_t obj(cache->line_cache_start);
      int sl,nr,nc,nra,nca,nrb,ncb,nrc,ncc,incc,incr,
        inca,incb;
      CBLAS_TRANSPOSE transa, transb;
      LineCache *lc, *lc0 = cache->line_cache_start;
      if (cache->nlines<2) goto NOT_SL;

      inca = lc0->inca;
      incb = lc0->incb;
      // Check all increments are equal
      for (int j=1; j<cache->nlines; j++) {
        lc = cache->line_cache_start+j;
        if (lc->inca != inca) goto NOT_SL;
        if (lc->incb != incb) goto NOT_SL;
      }

      // Check addresses in target (matrix result C) are RMO
      sl = obj.has_rmo(cache->nlines,superl_helper_t::c,
                  nrc,ncc,incr,incc);
      // rows in C must be contiguous
      if (!sl || incc!=1) goto NOT_SL;
      ldc = incr;

      // Check addresses in A are RMO
      sl = obj.has_rmo(cache->nlines,superl_helper_t::a,
                  nr,nc,incr,incc);
      if (!sl) goto NOT_SL;
      // Check whether A is transpose or not
      if (inca==1 && incc==0) {
        transa = CblasNoTrans;
        lda = incr; nra = nr; nca = nc;
      } else if (incr==0 && incc==1) {
        transa = CblasTrans;
        lda = inca; nra = nc; nca = nr;
      } else goto NOT_SL;

      // Check addresses in B are RMO
      sl = obj.has_rmo(cache->nlines,superl_helper_t::b,
                  nr,nc,incr,incc);
      if (!sl) goto NOT_SL;
      // Check whether B is transpose or not
      if (incb==1 && incc==0) {
        transb = CblasTrans;
        ldb = incr;
        nrb = nc; ncb = nr;
      } else if (incr==0 && incc==1) {
        transa = CblasNoTrans;
        ldb = incb;
        nrb = nr; ncb = nc;
      } else goto NOT_SL;

      // Verify matrix dimensions are OK
      assert(nca==nrb);
      assert(nra==nrc);
      assert(ncb==ncc);

      printf("check B: sl %d, nr %d, nc %d, incc %d, incr %d\n",
             sl,nr,nc,incc,incr);

      printf("check C: sl %d, nr %d, nc %d, incc %d, incr %d\n",
             sl,nr,nc,incc,incr);

      printf("check A: sl %d, nr %d, nc %d, incc %d, incr %d\n",
             sl,nr,nc,incc,incr);

      superlinear=1;
#define SHTRANS(trans) (trans==CblasTrans? 1 : 0)
      printf("superlinear with:\n"
             "A: transp %d, nr %d, nc %d, ld %d\n"
             "B: transp %d, nr %d, nc %d, ld %d\n"
             "C: ld %d\n",
             SHTRANS(transa),nra,nca,lda,
             SHTRANS(transb),nrb,ncb,ldb,
             ldc);
    NOT_SL: 

      exit(0);

#if 0
      // Detect if operation is superlinear
      LineCache *lc0, *lc1, *lc;
      int inca, incb, incc;
      int &nlines = cache->nlines;
      if (nlines>1) {
        lc0 = cache->line_cache_start;
        paa0 = *lc0->starta;
        pbb0 = *lc0->startb;
        pcc0 = lc0->target;
        inca = lc0->inca;
        ldb = lc0->incb;
        // Try to determine leading size of matrices A and B
        for (int j=0; j<nlines; j++) {
          lc = cache->line_cache_start+j;
          if (*lc->starta != paa0) {
            lda = *lc->starta-paa0;
            ldc = lc->target - pcc0;
            na = j;
            break;
          }
        }
        if (lda<=0) goto NOT_SL;
        lc1 = cache->line_cache_start+1;
        incb = *lc1->startb - pbb0;
        incc = lc1->target - lc0->target;
        if (incc != 1) goto NOT_SL;
        
        if (nlines % na != 0) goto NOT_SL;
        nb = nlines/na;
        // Check that is truly superlinear
        int l=0;
        for (int j=0; j<na; j++) {
          for (int k=0; k<nb; k++) {
            lc = cache->line_cache_start+l;
            if (lc->inca != inca) goto NOT_SL;
            if (lc->incb != ldb) goto NOT_SL;
            if (*lc->starta != paa0 + j*lda) goto NOT_SL;
            if (*lc->startb != pbb0 + k*incb) goto NOT_SL;
            if (lc->target != pcc0 + j*ldc + k*incc) goto NOT_SL;
            l++;
          }
        }
        // Apparently `dgemm' only works with contiguous matrices
        if (inca!=1 || incb!=1) goto NOT_SL;
      }
      superlinear=1;
    NOT_SL: 
//       printf("superlinear with na %d, lda %d, nb %d, ldb %d, ldc %d\n",
//              na,lda,nb,ldb,ldc);
#endif
    
      psc = new prod_subcache;
      assert(psc);
      assert(!cache->sc);
      cache->sc = psc;
      psc->superlinear = superlinear;
      if (superlinear) {
        psc->nra = nra;
        psc->nca = nca;
        psc->lda = lda;
        psc->ncb = ncb;
        psc->ldb = ldb;
        psc->ldc = ldc;
        psc->paa0 = paa0;
        psc->pbb0 = pbb0;
        psc->pcc0 = pcc0;
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
