//<require "prep.pl"//>//
//<=$warn_dont_modify //>

//__INSERT_LICENSE__
//$Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$
#include <math.h>
#include <stdio.h>
#ifdef USE_MKL
#include <mkl_cblas.h>
#endif

#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fastlib2.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int mem_size(const Indx & indx) {
  int size=1;
  for (int j=0; j<indx.size(); j++) {
    size *= indx[j];
  }
  return size;
}

//<$genone=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::__NAME__(const FastMat2 & A __OTHER_ARGS__) {

  CTX2_CHECK("__NAME__",this,&A);

  if (!ctx->was_cached  ) {
    assert(A.defined);
    Indx Afdims,fdims;
    A.get_dims(Afdims);
    // int Afd = Afdims.size();
  
    if (!defined) {
      create_from_indx(Afdims);
    }
    get_dims(fdims);
    if (Afdims != fdims) {
      Afdims.print("Dims of in matrix:");
      fdims.print("Dims of out matrix: ");
      assert(0);
    }

    int ndims = dims.size();
    Indx indx(ndims,1);
    int flag=1;
    while (flag) {
      cache->to_elems.push_back(location(indx));
      cache->from_elems.push_back(A.location(indx));
      flag=inc(indx,fdims);
    } 
    cache->nelems = cache->to_elems.size();
    cache->pto = &*cache->to_elems.begin();
    cache->pfrom = &*cache->from_elems.begin();

    ctx->op_count.put += cache->nelems;
    ctx->op_count.get += cache->nelems;
    __COUNT_OPER__;
  }

  double **pto,**pfrom,**pto_end;
  pto = cache->pto;
  pto_end = pto + cache->nelems;
  pfrom = cache->pfrom;
  while (pto < pto_end) {
    __ELEM_OPERATIONS__;
  }
  if (!ctx->use_cache) delete cache;
  return *this;
}  
//EOF
_//>

//<# This line expand the previous template to several routines //>//
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// ALL THESE ROUTINES GENERATED FROM THE $genone TEMPLATE
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
//< genone_all()//>


//<$gensetel=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::__NAME__(const double val, INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("__NAME__",this);
    Indx indx;
    READ_ARG_LIST(arg,indx,INT_ARG_LIST_DEFAULT_VAL,EXIT2)
    ctx->check(indx);
  }
#endif
  FastMatCache *cache = ctx->step();

  if (!ctx->was_cached  ) {
    Indx indx,fdims;
    get_dims(fdims);
    int ndims = fdims.size();
    assert(ndims>0);
#ifdef USE_VAR_ARGS
    va_list ap;
    va_start(ap,val);
    read_int_list(ndims,ap,&indx);
#else
    READ_INT_ARG_LIST(indx);
    assert(ndims==indx.size());
#endif
    cache->to = location(indx);
    ctx->op_count.put += 1;
    __COUNT_OPER__;
  }

  __ELEM_OPERATIONS__;

  if (!ctx->use_cache) delete cache;
  return *this;

}  
//EOF
_//>

//< gen_setel_all(); //>

//<$SET_ARRAY_KERNEL=<<'//EOF';
  double **pto,**pto_end;
  pto = cache->pto;
  pto_end = pto + cache->nelems;
  while (pto < pto_end) {
    **pto++ = *from++;
  }
  if (!ctx->use_cache) delete cache;
  return *this;
//EOF
_//>

//<$set_array=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::set(const double *a) {

  CTX2_CHECK("set",this);

  if (!ctx->was_cached  ) {
    if (!defined) {
      cache->nelems = 0;
      return *this;
    }
    // assert(defined);
    Indx fdims;
    get_dims(fdims);
    int ndims = dims.size();
    Indx indx(ndims,1);
    int size = mem_size(fdims);
    cache->to_elems.resize(size);
    cache->nelems = size;
    int flag=1,j=0;
    while (flag) {
      cache->to_elems[j++] = location(indx);
      flag=inc(indx,fdims);
    } 
    cache->pto = &*cache->to_elems.begin();
    ctx->op_count.get += cache->nelems;
    ctx->op_count.get += cache->nelems;
  }
  const double *from = a;

  __SET_ARRAY_KERNEL__;

}  
//EOF
_//>

//< print template_subst($set_array); //>//

//<$set_from_newmat = <<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::set(const Matrix & A) {

  // This should require that the arguments
  // should be generic ptrs. since A is a Newmat Matrix, not a
  // FastMat2. 
  // CTX2_CHECK("__NAME__",this,&A);
  CTX2_CHECK("set",this);

  if (!ctx->was_cached  ) {
    int m = A.Nrows();
    int n = A.Ncols();

    Indx newdims,fdims;
    newdims.push_back(m);
    newdims.push_back(n);
  
    if (!defined) 
      create_from_indx(newdims);
    if (m==0 || n==0) {
      cache->nelems=0;
      return *this;
    }

    get_dims(fdims);
    assert(fdims == newdims);

    int size = mem_size(fdims);
    cache->nelems = size;
    cache->to_elems.resize(size);
    cache->pto = &*cache->to_elems.begin();

    int jj=0;
    Indx indx(2,0);
    for (int j=1; j<=m; j++) {
      indx[0]=j;
      for (int k=1; k<=n; k++) {
	indx[1]=k;
	cache->to_elems[jj++] = location(indx);
      }
    }
    // cache->from = A.Store();
  }
  double *from = A.Store();

  __SET_ARRAY_KERNEL__;

}
//EOF
_//>

//< print template_subst($set_from_newmat); //>//

//<$DEFINE_SUB_INDICES_FROM_ARGS=<<'//EOF';
    sindx.push_back(m);
    va_list ap;
    va_start(ap,m);
    read_int_list(Afdims.size()-1,ap,&sindx);

    int nfree=0,nc=0;
    for (int j=0; j<sindx.size(); j++) {
      int k = sindx[j];
      if (k>0 && k>nfree) nfree = k;
      if (k<0) nc++;
    }

    Indx ifree(nfree,0),icontr(nc,0);
    int ic=0;
    for (int j=0; j<sindx.size(); j++) {
      int k = sindx[j];
      if (k>0) {
	ifree[k-1] = j+1;
      } else {
	icontr[ic++] = j+1;
      }
    }
//EOF
_//>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
//<$DEFINE_SUB_INDICES_ALL=<<'//EOF';
sindx = Indx(Afdims.size(),-1);

int nfree=0,nc=0;
for (int j=0; j<sindx.size(); j++) {
  int k = sindx[j];
  if (k>0 && k>nfree) nfree = k;
  if (k<0) nc++;
}

Indx ifree(nfree,0),icontr(nc,0);
int ic=0;
for (int j=0; j<sindx.size(); j++) {
  int k = sindx[j];
  if (k>0) {
    ifree[k-1] = j+1;
  } else {
    icontr[ic++] = j+1;
  }
}
//EOF
_//>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// This cache is used when converting a "generic-sum"
// function like `norm_p' or `sum' to the corresponding
// `_all' version, for instance `norm_p' -> `norm_p_all'. 
// In that case we need a sub_cache that stores a FastMat2 
// that is a scalar. 
class gensum_all_cache : public FastMatSubCache {
public:
  FastMat2 tmp;
  gensum_all_cache(FastMat2::CacheCtx *ctxp) 
    : tmp(ctxp) { }
  ~gensum_all_cache() {};
};

//<$gen_sum=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/* Obtained from pattern $gen_sum with args;
   'INI_LOOP' => '__INI_LOOP__'
   'NAME' => '__NAME__'
   'ELEM_OPERATIONS' => '__ELEM_OPERATIONS__'
   'COUNT_OPER' => '__COUNT_OPER__'
   'OTHER_ARGS' => '__OTHER_ARGS__'
   'C' => '__C__'
   'POST_LOOP_OPS' => '__POST_LOOP_OPS__'
*/
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::__NAME__(const FastMat2 & A, __OTHER_ARGS__ __C__ 
			      const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("__NAME__",this,&A);
    Indx indx;
    indx.push_back(m);
    READ_ARG_LIST(arg,indx,INT_ARG_LIST_DEFAULT_VAL,EXIT2)
    ctx->check(indx);
  }
#endif
  FastMatCache *cache = ctx->step();

  if (!ctx->was_cached  ) {
    Indx sindx,fdims,Afdims;
    assert(A.defined);
    A.get_dims(Afdims);

    if (m!=0) {
      sindx.push_back(m);
#ifdef USE_VAR_ARGS
      va_list ap;
      va_start(ap,m);
      read_int_list(Afdims.size()-1,ap,&sindx);
#else
      READ_INT_ARG_LIST(sindx);
      assert(sindx.size() == Afdims.size());
#endif
    } else {
      sindx = Indx(Afdims.size(),-1);
    }

    int nfree=0,nc=0;
    for (int j=0; j<sindx.size(); j++) {
      int k = sindx[j];
      if (k>0 && k>nfree) nfree = k;
      if (k<0) nc++;
    }

    Indx ifree(nfree,0),icontr(nc,0);
    int ic=0;
    for (int j=0; j<sindx.size(); j++) {
      int k = sindx[j];
      if (k>0) {
	ifree[k-1] = j+1;
      } else {
	icontr[ic++] = j+1;
      }
    }
  
    Indx ndimsf(nfree,0),ndimsc(nc,0);
    int nlines=1;
    for (int j=0; j<nfree; j++) {
      int k = ifree[j];
      ndimsf[j] = Afdims[k-1];
      nlines *= ndimsf[j];
    }

    // Dimension B (*this) if necessary
    if (!defined) {
      create_from_indx(ndimsf);
    }

    get_dims(fdims);
    assert(ndimsf == fdims);

    int line_size=1;
    for (int j=0; j<nc; j++) {
      int k = icontr[j];
      ndimsc[j] = Afdims[k-1];
      line_size *= ndimsc[j];
    }

    int ndims = Afdims.size();
    Indx findx(nfree,1),cindx(nc,1),tot_indx(ndims,0);

    // Loading addresses in cache
    // For each element in the distination target, we store the complete
    // list of addresses of the lines of elements that contribute to
    // it. 
    cache->prod_cache.resize(nlines);
    cache->line_cache_start = &*cache->prod_cache.begin();
    cache->nlines = nlines;
    cache->line_size = line_size;
    LineCache *lc;
    for (int jlc=0; jlc<nlines; jlc++) {
      lc = cache->line_cache_start + jlc;
      lc->linea.resize(line_size);
      // cache->prod_cache.push_back(LineCache());
      lc->target = location(findx);
      // findx.print("for free indx: ");

      cindx= Indx(nc,1);
      for (int j=0; j<nfree; j++)
	tot_indx[ifree[j]-1] = findx[j];

      int kk=0;
      while(1) {
	for (int j=0; j<nc; j++) {
	  int k=icontr[j];
	  tot_indx[k-1] = cindx[j];
	}
	// tot_indx.print("tot_indx: ");

	lc->linea[kk++] = A.location(tot_indx);
	if (!inc(cindx,ndimsc)) break;
      }
      lc->starta = &*lc->linea.begin();
      if (!inc(findx,ndimsf)) break;
    }
    int ntot = nlines*line_size;
    ctx->op_count.get += ntot;
    ctx->op_count.put += nlines;
    __COUNT_OPER__;
  }

  LineCache *lc;
  double **pa,**pe,val;
  __PRE_LOOP_OPS__;
  for (int j=0; j<cache->nlines; j++) {
    lc = cache->line_cache_start+j;
    pa = lc->starta;
    pe = pa + cache->line_size;
    // val=0;
    __INI_LOOP__;
    while (pa<pe) {
      __ELEM_OPERATIONS__;
    }
    __POST_LOOP_OPS__;
    *lc->target = val;
  }
  __AFTER_ALL_STRIDES__;
  if (!ctx->use_cache) delete cache;
  return *this;
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double FastMat2::__NAME___all(__OTHER_ARGS__) const {

  CTX2_CHECK("__NAME___all",this);

  gensum_all_cache *gsac=NULL;
  if (!ctx->was_cached) {
    gsac = new gensum_all_cache(ctx);
    assert(gsac);
    assert(!cache->sc);
    cache->sc = gsac;
  }

  gsac = dynamic_cast<gensum_all_cache *> (cache->sc);
  assert(gsac);

  double retval;
  FastMat2 &tmp = gsac->tmp;
  tmp.__NAME__(*this __C__ __OTHER_ARGS_P__);
  retval = double(tmp);
  if (!ctx->use_cache) delete cache;
  return retval;
}

//EOF
_//>
//< gen_sum_all(); //>

//<$in_place = <<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::__NAME__(__FUN_ARGS__) {

  CTX2_CHECK("__NAME__",this);

  if (!ctx->was_cached  ) {
    assert(defined);

    Indx fdims;
    get_dims(fdims);

    int ndims = dims.size();
    Indx indx(ndims,1);
    int flag=1;
    cache->nelems = mem_size(fdims);
    cache->to_elems.resize(cache->nelems);
    int j=0;
    while (flag) {
      cache->to_elems[j++] = location(indx);
      flag=inc(indx,fdims);
    } 
    cache->pto = &*cache->to_elems.begin();
  }

  double **to = cache->pto;
  double **to_end = to + cache->nelems;
  while (to<to_end) {
    __ELEM_OPERATIONS__;
  }

  if (!ctx->use_cache) delete cache;
  return *this;
}  
//EOF
_//>
//< in_place_all(); //>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

//<$DEFINE_NEWMAT=<<'//EOF';
assert(ndims<=2);
int m=(ndims>=1 ?  dims[0].dim : 1);
int n=(ndims==2 ?  dims[1].dim : 1);
if (A.Nrows()==0) {
  A.ReSize(m,n);
}
//EOF
_//>


//<$export_array=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
__CONST__ FastMat2 & FastMat2::export_vals(__ARG__) __CONST__ {

  CTX2_CHECK("export_vals",this);

  if (!ctx->was_cached  ) {
    if (!defined) {
      cache->nelems = 0;
      return *this;
    }
    // assert(defined);
    Indx fdims;
    get_dims(fdims);
    int ndims=fdims.size();

    __DEFINE_NEWMAT_MAYBE__;

    Indx indx(ndims,1);
    int size = mem_size(fdims);
    cache->from_elems.resize(size);
    cache->nelems = size;
    int flag=1,j=0;
    while (flag) {
      cache->from_elems[j++] = location(indx);
      flag=inc(indx,fdims);
    } 
    cache->pfrom = &*cache->from_elems.begin();
  }
  __DEFINE_TARGET_POINTER__;

  double **pfrom,**pfrom_end;
  pfrom = cache->pfrom;
  pfrom_end = pfrom + cache->nelems;
  while (pfrom < pfrom_end) {
    *to++ = **pfrom++;
  }
  if (!ctx->use_cache) delete cache;
  return *this;
}  
//EOF
_//>

//< export_vals_array();//>

//<$contraction=<<'//EOF';
FastMat2 & FastMat2::ctr(const FastMat2 & A,
                         const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("__NAME__",this,&A);
    ctx->check(m);
    Indx indx;
    READ_ARG_LIST(arg,indx,INT_ARG_LIST_DEFAULT_VAL,EXIT2)
    ctx->check(indx);
  }
#endif
  FastMatCache *cache = ctx->step();

  if (!ctx->was_cached  ) {
    Indx ia,ii;
  
    Indx Afdims,fdims;
    assert(A.defined);
    A.get_dims(Afdims);

    // maxc:= maximum contracted index
    ii.push_back(m);
    int ndims = Afdims.size();

#ifdef USE_VAR_ARGS
    va_list ap;
    va_start(ap,m);
    read_int_list(ndims-1,ap,&ii);
    va_end(ap);
#else
    READ_INT_ARG_LIST(ii);
    assert(ii.size() == ndims);
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
  
    Indx ndimsf(nfree,0),ndimsc(nc,0);
    for (int j=0; j<nfree; j++) {
      int k = ifree[j];
      ndimsf[j] = Afdims[k-1];
    }

    // Dimension B (*this) if necessary
    if (!defined) {
      create_from_indx(ndimsf);
    }

    get_dims(fdims);
    assert(ndimsf == fdims);

    for (int j=0; j<nc; j++) {
      int k1 = icontr[2*j];
      int k2 = icontr[2*j+1];
      int PFUNUSED nd1, nd2;
      nd1 = Afdims[k1-1];
      nd2 = Afdims[k2-1];
      assert(nd1==nd2);
      ndimsc[j]=nd1;
    }

    // ndimsc.print("dimensions of the contracted part: ");

    Indx findx(nfree,1),cindx(nc,1),tot_indx(ndims,0);

    // Loading addresses in cache
    // For each element in the distination target, we store the complete
    // list of addresses of the lines of elements that contribute to
    // it. 
    int nlines = mem_size(ndimsf);
    int line_size = mem_size(ndimsc);
    cache->prod_cache.resize(nlines);
    cache->nlines = nlines;
    cache->line_size = line_size;
    cache->line_cache_start = &*cache->prod_cache.begin();

    int jlc=0;
    while (1) {
      LineCache *lc = cache->line_cache_start + jlc++;
      lc->linea.resize(line_size);
      lc->starta = &*lc->linea.begin();
      lc->target = location(findx);

      cindx= Indx(nc,1);
      for (int j=0; j<nfree; j++)
	tot_indx[ifree[j]-1] = findx[j];

      int jj=0;
      while(1) {
      
	for (int j=0; j<nc; j++) {
	  int k1=icontr[2*j];
	  int k2=icontr[2*j+1];
	  tot_indx[k1-1] = cindx[j];
	  tot_indx[k2-1] = cindx[j];
	}
	// tot_indx.print("tot_indx: ");

	lc->linea[jj++] = A.location(tot_indx);
	if (!inc(cindx,ndimsc)) break;
      }
      if (!inc(findx,ndimsf)) break;
    }
    int ntot = cache->nlines*cache->line_size;
    ctx->op_count.get += ntot;
    ctx->op_count.put += cache->nlines;
    ctx->op_count.sum += ntot;
  }

  // Perform computations using cached addresses
  int nlines = cache->nlines;
  double **pa,**pa_end,sum;
  for (int j=0; j<nlines; j++) {
    LineCache *lc = cache->line_cache_start+j;
    pa = lc->starta;
    pa_end = pa + cache->line_size;
    sum=0.;
    while (pa<pa_end) {
      sum += **pa++;
    }
    *(lc->target) = sum;
  }

  if (!ctx->use_cache) delete cache;
  return *this;
}
//EOF
_//>

//< print template_subst($contraction); //>//

//<$diag=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::diag(FastMat2 & A,const int m,
                          INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("__NAME__",this,&A);
    ctx->check(m);
    Indx indx;
    READ_ARG_LIST(arg,indx,INT_ARG_LIST_DEFAULT_VAL,EXIT2)
    ctx->check(indx);
  }
#endif
  FastMatCache *cache = ctx->step();

  if (!ctx->was_cached) {
    Indx ia,ii;
  
    Indx Afdims;
    assert(A.defined);
    A.get_dims(Afdims);

    ii.push_back(m);
    int ndims = Afdims.size();

#ifdef USE_VAR_ARGS
    va_list ap;
    va_start(ap,m);
    read_int_list(ndims-1,ap,&ii);
    va_end(ap);
#else
    READ_INT_ARG_LIST(ii);
    assert(ndims==ii.size());
#endif

    int nfree=0,nc=0;
    Indx indxmap = ii;
    for (int j=0; j<ii.size(); j++) {
      int k = ii[j];
      if (k>0) nfree++;
      if (k<0) {
	nc++;
	indxmap[j]=-k;
      }
    }
    nc /= 2;
    int ntot = nfree+nc;

    Indx newdims(ntot,0);
    for (int j=0; j<ndims; j++) {
      int k = indxmap[j];
      newdims[k-1] = Afdims[indxmap[j]-1];
    }

    // Dimension B (*this) if necessary
    if (!defined) {
      create_from_indx(newdims);
    }

    Indx fdims;
    get_dims(fdims);
    assert(newdims == fdims);

  // B.diag(A,indices...)
    Indx bindx(ntot,1),aindx(ndims,0);

    // Loading addresses in cache
    // For each element in the distination target, we store the complete
    // list of addresses of the lines of elements that contribute to
    // it. 
    int nelems = mem_size(newdims);
    cache->to_elems.resize(nelems);
    cache->from_elems.resize(nelems);
    cache->nelems = nelems;
    cache->pto = &*cache->to_elems.begin();
    cache->pfrom = &*cache->from_elems.begin();
    int jj=0;
    while (1) {
      for (int k=0; k<ndims; k++) 
	aindx[k] = bindx[indxmap[k]-1];
    
      cache->to_elems[jj] = location(bindx);
      cache->from_elems[jj] = A.location(aindx);
      jj++;
      if (!inc(bindx,fdims)) break;
    }
    ctx->op_count.get += nelems;
    ctx->op_count.put += nelems;
  }

  // Perform computations using cached addresses
  double **pto, **pfrom, **pto_end;
  pto = cache->pto;
  pto_end = cache->pto + cache->nelems;
  pfrom = cache->pfrom;
  
  while  (pto<pto_end) {
    **pto++ = **pfrom++;
  }

  if (!ctx->use_cache) delete cache;
  return *this;
}
//EOF
_//>

//< print template_subst($diag); //>//

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

//<$get=<<'//EOF';
double FastMat2::get(INT_VAR_ARGS_ND) const {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("__NAME__",this);
    Indx indx;
    READ_ARG_LIST(arg,indx,INT_ARG_LIST_DEFAULT_VAL,EXIT2)
    ctx->check(indx);
  }
#endif
  FastMatCache *cache = ctx->step();

  if (!ctx->was_cached) {
    Indx indx,fdims;
    get_dims(fdims);
    int ndims = fdims.size();

#ifdef USE_VAR_ARGS
    assert(0); // No funciona mas con variable args
    va_list ap;
    va_start(ap,i);
    read_int_list(ndims-1,ap,&indx);
    va_end(ap);
#else
    READ_INT_ARG_LIST(indx);
    assert(indx.size()==ndims);
#endif

    cache->from = location(indx);
    ctx->op_count.get += 1;
  }

  if (!ctx->use_cache) delete cache;
  return *cache->from;
}
//EOF
_//>

//< print template_subst($get); //>//

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

//<$operator_double=<<'//EOF';
// This is somewhat redundant because you can use A.get()
FastMat2::operator double() const { 
  
  CTX2_CHECK("double",this);

  if (!ctx->was_cached) {
    Indx fdims;
    get_dims(fdims);
    assert(fdims.size()==0);
    ctx->op_count.get += 1;
    ctx->op_count.put += 1;
  }

  if (!ctx->use_cache) delete cache;
  return *store;
}
//EOF
_//>

//< print template_subst($operator_double); //>//

//<$det=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double FastMat2::det(void) const {

  CTX2_CHECK("det",this);

  if (!ctx->was_cached) {
    Indx dims_;
    get_dims(dims_);
    int ndims = dims_.size();
    assert (ndims==2);
    assert (dims_[0] == dims_[1]);
    int m= dims_[0];
    if (m<=3) {
      cache->nelems=m*m;
      cache->from_elems.resize(cache->nelems);
      cache->pfrom = &*cache->from_elems.begin();
      int jj=0;
      Indx indx(2,0);
      for (int j=1; j<=m; j++) {
	for (int k=1; k<=m; k++) {
	  indx[0]=j;
	  indx[1]=k;
	  cache->pfrom[jj++] = location(indx);
	}
      }
    } else if (m>3) {
      cache->A = new Matrix(m,m);
    }
  }

  int m = this->dims_p[0].dim;
  double **pfrom = cache->pfrom, det_;
#define A(i,j) (*pfrom[(i-1)*M+(j-1)])
  if (m==1) {
    det_ = **pfrom;
  } else if(m==2) {
#define M 2
    det_ = A(1,1)*A(2,2)-A(1,2)*A(2,1);
#undef M
  } else if (m==3) {
#define M 3
    det_ = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      + A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))
      + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1));
#undef M
#undef A
  } else {
    export_vals(*cache->A);
    LogAndSign ld;
    ld = cache->A->LogDeterminant();
    det_ = ld.Value();
  }
  if (!ctx->use_cache) delete cache;
  return det_;
}
//EOF
_//>

//< print template_subst($det); //>//

//<$kron=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::kron(const FastMat2 & A,
                          const FastMat2 & B) {
  CTX2_CHECK("det",this,&A,&B);

  if (!ctx->was_cached) {
    Indx Adims,Bdims,dims_;
    A.get_dims(Adims);
    B.get_dims(Bdims);
    int ndims = Adims.size();
    assert(Bdims.size()==ndims);
    
    Indx newdims(ndims,0);
    for (int j=0; j<ndims; j++) 
      newdims[j] = Adims[j]*Bdims[j];

    if(!defined)
      create_from_indx(newdims);

    get_dims(dims_);
    assert(dims_==newdims);

    int nelems = mem_size(newdims);
    // Here we have nelems elements in A B and *this.
    // We store the A and B addresses in from_elems and the *this
    // adresses in to_elems 
    cache->nelems = nelems;
    cache->from_elems.resize(2*nelems);
    cache->to_elems.resize(nelems);
    cache->pto = &*cache->to_elems.begin();
    cache->pfrom = &*cache->from_elems.begin();
    double **pfrom_b = cache->pfrom + nelems;
    double **pfrom_a = cache->pfrom;
    double **pto = cache->pto;

    // Load addresses in cache
    Indx Aindx(ndims,1),Bindx(ndims,1),indx_(ndims,0);
    int jj=0;
    while (1) {
      while (1) {
	pfrom_a[jj] = A.location(Aindx);
	pfrom_b[jj] = B.location(Bindx);

	for (int jd=0; jd<ndims; jd++) 
	  indx_[jd] = (Aindx[jd]-1)*Bdims[jd] + Bindx[jd];
	pto[jj] = location(indx_);

	jj++;
	if (!inc(Bindx,Bdims)) break;
      }
      if (!inc(Aindx,Adims)) break;
    }

    ctx->op_count.get += nelems;
    ctx->op_count.put += nelems;
    ctx->op_count.mult += nelems;
  }

  // Perform computations using cached addresses
  double **pto, **pfrom_a, **pfrom_b, **pto_end;
  pto = cache->pto;
  pto_end = cache->pto + cache->nelems;
  pfrom_a = cache->pfrom;
  pfrom_b = cache->pfrom + cache->nelems;
  
  while  (pto<pto_end) 
    (**pto++) = (**pfrom_a++) * (**pfrom_b++);

  if (!ctx->use_cache) delete cache;
  return *this;

}    
//EOF
_//>

//< print template_subst($kron); //>//


//<$eye=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::eye(const double a) {
  set(0.);

  CTX2_CHECK("eye",this);

  if (!ctx->was_cached  ) {

    assert(defined);

    Indx fdims;
    get_dims(fdims);
    int ndims = fdims.size();
    assert(ndims==2);
    assert(fdims[0]==fdims[1]);
    int n = fdims[0];

    Indx indx;
    cache->to_elems.resize(n);
    for (int j=1; j<=n; j++) {
      indx = Indx(2,j);
      cache->to_elems[j-1] = location(indx);
    } 
    cache->nelems = n;

    cache->pto = &*cache->to_elems.begin();

    ctx->op_count.put += n;
    __COUNT_OPER__;
  }

  double **pto = cache->pto;
  double **pto_end = pto + cache->nelems;

  while (pto < pto_end) {
    **pto++ = a;
  }

  if (!ctx->use_cache) delete cache;
  return *this;
}  
//EOF
_//>

//< print template_subst($eye); //>//

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class detsur_cache : public FastMatSubCache {
public:
  FastMat2 g;
  int m,n;
  double *v[6], *nor_p[3];
  detsur_cache(FastMat2::CacheCtx *ctxp) 
    : g(ctxp) { }
  ~detsur_cache() { }
};

//<$detsur=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double FastMat2::detsur(FastMat2 *nor) {

  CTX2_CHECK("detsur",this,nor);

  detsur_cache * dsc;
  if (!ctx->was_cached) {
    Indx fdims;
    get_dims(fdims);
    assert(fdims.size()==2);

    dsc = new detsur_cache(ctx);
    assert(dsc);
    dsc->m = dim(1);
    dsc->n = dim(2);
    if (dsc->m>0) dsc->g.resize(2,dsc->m,dsc->m);
    assert(!cache->sc);
    cache->sc = dsc;
    if (nor) {
      // The jacobian should be n-1 x n
      assert(dsc->m == dsc->n-1);
      // only dimensions 2 and 3 are allowed
      assert(dsc->n==2 || dsc->n==3);
      if (nor->defined) {
	// Check the normal has to be a vector
	assert(nor->n()==1);
	// the dimension of nor has to be of size #n#
	assert(nor->dim(1)==dsc->n);
      } else {
	// Create a vector of length dsc->n
	Indx nor_dims;
	nor_dims.push_back(dsc->n);
	nor->create_from_indx(nor_dims);
      }
      // cached vector of pointers to elements in the Jacobian
      double **w = dsc->v;
      Indx indx,indx_n;
      for (int k=1; k<=dsc->n-1; k++) {
	for (int l=1; l<=dsc->n; l++) {
	  indx.reset();
	  indx.push_back(k);
	  indx.push_back(l);
	  *w++ = location(indx);
	}
      }

      // cached vector of pointers to elements in normal vector
      w = dsc->nor_p;
      for (int l=1; l<=dsc->n; l++) {
	indx.reset();
	indx.push_back(l);
	*w++ = nor->location(indx);
      }
    }
  }

  dsc = dynamic_cast<detsur_cache *> (cache->sc);
  double retval;
  if (dsc->m == 0) return 1.;
  if (nor) {
#define VV(j) (*vv[j])
#define NN(j) (*nn[j])
    double **nn = dsc->nor_p;
    double **vv = dsc->v;
    if (dsc->n==2) {
      *nn[0] = +VV(1);
      *nn[1] = -VV(0);
      retval = sqrt(NN(0)*NN(0)+NN(1)*NN(1));
    } else {
#define JACO(j,k) (*vv[((j)-1)*3+(k)-1])
      NN(0) = -JACO(1,2)*JACO(2,3)+JACO(1,3)*JACO(2,2);
      NN(1) = -JACO(1,3)*JACO(2,1)+JACO(1,1)*JACO(2,3);
      NN(2) = -JACO(1,1)*JACO(2,2)+JACO(1,2)*JACO(2,1);
      retval = sqrt(NN(0)*NN(0)+NN(1)*NN(1)+NN(2)*NN(2));
    }
  } else {
    dsc->g.prod(*this,*this,1,-1,2,-1);
    retval = sqrt(dsc->g.det());
  }
  if (!ctx->use_cache) delete cache;
  return retval;
}
//EOF
_//>

//< print template_subst($detsur); //>//

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class cross_cache : public FastMatSubCache {
public:
  int ndim;
  const double **a,**b;
  double **c;
  cross_cache() { a=b=NULL; c=NULL; }
  ~cross_cache();
};

//<$cross=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::cross(const FastMat2 & a,const FastMat2 & b) {

  CTX2_CHECK("cross",this,&a,&b);

  // Cross product of vectors
  cross_cache *ccache;
  if (!ctx->was_cached) {
    assert(a.defined);
    assert(a.n()==1);
    assert(b.defined);
    assert(b.n()==1);
    assert(a.dim(1)==b.dim(1));
    assert(a.dim(1)==3 || a.dim(1)==2);
    // Number of space dimensions
    int ndim = a.dim(1);
    // Length of out vector (Cross Product Size)
    // In 2D it is a scalar
    int cps = (ndim==3 ? 3 : 1);
    if (defined) {
      assert(n()==1);
      assert(dim(1)==cps);
    } else resize(1,cps);
    
    ccache = new cross_cache();
    assert(ccache);
    ccache->ndim = ndim;
    ccache->a = new const double *[ndim];
    ccache->b = new const double *[ndim];
    ccache->c = new double *[cps];

    Indx indx(1,0);
    for (int j=1; j<=ndim; j++) {
      indx[0]=j;
      ccache->a[j-1] = a.location(indx);
      ccache->b[j-1] = b.location(indx);
    }
    indx = Indx(1,0);
    for (int j=1; j<=cps; j++) {
      indx[0]=j;
      ccache->c[j-1] = location(indx);
    }
    assert(!cache->sc);
    cache->sc = ccache;
  }

  ccache  = dynamic_cast<cross_cache *> (cache->sc);
  const double **&aa = ccache->a;
  const double **&bb = ccache->b;
  double **&c = ccache->c;

  if (ccache->ndim==3) {
    *c[0] = *aa[1] * *bb[2] - *aa[2] * *bb[1];
    *c[1] = *aa[2] * *bb[0] - *aa[0] * *bb[2];
    *c[2] = *aa[0] * *bb[1] - *aa[1] * *bb[0];
  } else {
    *c[0] = *aa[0] * *bb[1] - *aa[1] * *bb[0];
  }    

  if (!ctx->use_cache) delete cache;
  return *this;
}
//EOF
_//>

//< print template_subst($cross); //>//

cross_cache::~cross_cache() {
  delete[] a;
  delete[] b;
  delete[] c;
}

//<=$warn_dont_modify //>
