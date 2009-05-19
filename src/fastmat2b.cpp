//__INSERT_LICENSE__
#include <cmath>
#include <ctime>
#include <cstdio>
#include <unistd.h>

#include <algorithm>
#include <mkl_cblas.h>

#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fm2stats.h>
#include <src/fastlib2.h>

FastMat2Stats glob_fm2stats;
int FASTMAT2_USE_DGEMM=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class prod_subcache_t : public FastMatSubCache {
public:
  int superlinear, lda, ldb, ldc, nra, nca, ncb,
    not_superlinear_ok_flag;
  CBLAS_TRANSPOSE transa, transb;
  double *paa0, *pbb0, *pcc0;
  enum  mode_t { a,b,c,none };
  FastMatCache *cache;
  LineCache *line_cache_start;
  prod_subcache_t(FastMatCache *cache_a)
    : superlinear(0),
      cache(cache_a), 
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
  void ident();
  void dgemm();
  void print();
  void ok() { not_superlinear_ok_flag=1; }
  int not_superlinear_ok();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void prod_subcache_t::ident() {

  not_superlinear_ok_flag=0;
  superlinear = 0;
  int N = cache->nlines;

  if (N<=2) { ok(); return; }
  double 
    *aa00 = address(0,prod_subcache_t::a),
    *bb00 = address(0,prod_subcache_t::b),
    *cc00 = address(0,prod_subcache_t::c);
  
  nra = -1;
  nca = cache->line_size;
  if (nca==1) { ok(); return; }
  int j,taflag=0,tbflag=0;

  for (j=1; j<N; j++) 
    if ((address(j,prod_subcache_t::a) - aa00)!=0) break;
  ncb = j;
  if (N % ncb != 0) return;
  nra = N/ncb;

  LineCache *lc0 = cache->line_cache_start;
  int incb,incc,inca = lc0->inca;
  ldb = lc0->incb;
  incc = address(1,prod_subcache_t::c) - cc00;

  if (nra>1) {
    lda = address(ncb,prod_subcache_t::a) - aa00,
    ldc = address(ncb,prod_subcache_t::c) - cc00;
  } else {
    lda = (inca!=1? 1 : nca);
    ldc = (incc!=1? 1 : ncb);
  }

  if (ncb>1) {
    incb = address(1,prod_subcache_t::b) - bb00,
    incc = address(1,prod_subcache_t::c) - cc00;
  } else {
    incb = 1;
    incc = 1;
  }

  int dp, l=0;
  for (j=0; j<nra; j++) {
    for (int k=0; k<ncb; k++) {
      if (lc0[l].inca != inca) return;
      if (lc0[l].incb != ldb) return;

      dp = address(l,prod_subcache_t::a) - aa00; 
      if (dp != lda*j) return;

      dp = address(l,prod_subcache_t::b) - bb00; 
      if (dp != incb*k) return;

      dp = address(l,prod_subcache_t::c) - cc00; 
      if (dp != incc*k + ldc*j) return;

      l++;
    }
  }

  if (incc!=1) return;

  if (inca!=1) {
    if (lda==1) {
      swap(lda,inca);
      taflag=1;
    } else return;
  }

  if (incb!=1) {
    if (ldb==1) {
      swap(ldb,incb);
      tbflag = 1;
    } else return;
  }

  if (lda<nca || ldb<ncb || ldc<ncb) return;

  transa = (taflag? CblasTrans : CblasNoTrans);
  transb = (tbflag? CblasTrans : CblasNoTrans);

  paa0 = *lc0->starta;
  pbb0 = *lc0->startb;
  pcc0 = lc0->target;

  superlinear = 1;
  return;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int prod_subcache_t::not_superlinear_ok() {
  return superlinear || not_superlinear_ok_flag || 
    nra==1 || nca==1 || ncb==1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void prod_subcache_t::print() {
  printf("dgemm args: sl %d,transa %d,transb %d,nra %d,nca %d,"
         "ncb %d,lda %d,ldb %d,ldc %d\n",
         superlinear,transa==CblasTrans,transb==CblasTrans,
         nra,nca,ncb,lda,ldb,ldc);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void prod_subcache_t::dgemm() {
#ifdef USE_MKL
  cblas_dgemm(CblasRowMajor,transa,transb,
              nra,ncb,nca,1.0,paa0,lda,pbb0,ldb,0.0,
              pcc0,ldc);
#else
  abort();
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & FastMat2::prod(const FastMat2 & A,const FastMat2 & B,
                          const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod",this,&A,&B);
    ctx->check(m);
    Indx indx;
    READ_ARG_LIST(arg,indx,INT_ARG_LIST_DEFAULT_VAL,EXIT2)
    ctx->check(indx);
  }
#endif
  FastMatCache *cache = ctx->step();

  prod_subcache_t *psc=NULL;
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

    psc = new prod_subcache_t(cache);
    assert(psc);
    assert(!cache->sc);
    cache->sc = psc;
    psc->ident();
    if (FASTMAT2_USE_DGEMM) {
      if (!psc->not_superlinear_ok()) {
        printf("NOT SL!!\n");
      }
      // psc->print();
      
      glob_fm2stats.was_sl_count += psc->superlinear;
      glob_fm2stats.was_not_sl_count += !psc->superlinear;
      
      if (psc->superlinear) {
        glob_fm2stats.was_sl= 1;
        psc->superlinear = glob_fm2stats.use_dgemm;
      }
    }
  }

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
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2Stats::report() {
  int total = 
    was_not_sl_count + was_sl_count;
  double ratio = 100.0*double(was_sl_count)/total;
//   printf("FM2STATS: sl %d(%.3f%%), not sl %d, total %d\n",
//          was_sl_count,ratio,was_not_sl_count,total);
  FILE *fid = fopen("/tmp/fm2stats.log","a");
  char *cwd = getcwd(NULL,0);
  time_t tt = time(NULL);
  // char *t = asctime(localtime(&tt));
#define MXTM 1000
  char t[MXTM];
  strftime(t,MXTM,"%a, %d %b %Y %H:%M:%S %z",
           localtime(&tt));
  fprintf(fid,"FM2STATS: [%s] sl %d(%.3f%%), not sl %d, "
          "total %d, cwd %s\n",
          t,was_sl_count,ratio,was_not_sl_count,total,cwd); 
  free(cwd);
  fclose(fid);
}
