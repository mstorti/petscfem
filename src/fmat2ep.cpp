//<require "prep.pl"//>//
//<=$warn_dont_modify //>

//__INSERT_LICENSE__
//$Id: fmat2ep.cpp,v 1.9 2001/06/08 14:25:56 mstorti Exp $
#include <math.h>
#include <stdio.h>

#include "fem.h"
#include "fastmat2.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int mem_size(const Indx & indx) {
  int size=1;
  for (int j=0; j<indx.size(); j++) {
    size *= indx[j];
  }
  return size;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

//<$cache_op=<<'//EOF';
FastMatCache *cache;
LineCache *lc;

if (was_cached) {
  cache = cache_list_begin[position_in_cache++];
#ifdef FM2_CACHE_DBG
  printf ("reusing cache: ");
#endif
} else if (!use_cache) {
  cache = new FastMatCache;
} else {
  cache = new FastMatCache;
  cache_list->push_back(cache);
  cache_list_begin = cache_list->begin();
  cache_list->list_size =
    cache_list_size = cache_list->size();
  position_in_cache++;
#ifdef FM2_CACHE_DBG
  printf ("defining cache: ");
#endif
}
#ifdef FM2_CACHE_DBG
printf(" cache_list %p, cache %p, position_in_cache %d\n",
       cache_list,cache,position_in_cache-1);
#endif
//EOF
_//>
//<$CACHE_OPERATIONS = $cache_op;//>

//<$genone=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::__NAME__(const FastMat2 & A __OTHER_ARGS__) {

  __CACHE_OPERATIONS__;

  if (!was_cached  ) {
    assert(A.defined);
    Indx Afdims,fdims;
    A.get_dims(Afdims);
    int Afd = Afdims.size();
  
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
    cache->pto = cache->to_elems.begin();
    cache->pfrom = cache->from_elems.begin();

    op_count.put += cache->nelems;
    op_count.get += cache->nelems;
    __COUNT_OPER__;
  }

  double **pto,**pfrom,**pto_end;
  pto = cache->pto;
  pto_end = pto + cache->nelems;
  pfrom = cache->pfrom;
  while (pto < pto_end) {
    __ELEM_OPERATIONS__;
  }
  if (!use_cache) delete cache;
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
FastMat2 & FastMat2::__NAME__(const double val, INT_VAR_ARGS) {

  __CACHE_OPERATIONS__;

  if (!was_cached  ) {
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
    op_count.put += 1;
    __COUNT_OPER__;
  }

  __ELEM_OPERATIONS__;

  if (!use_cache) delete cache;
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
  if (!use_cache) delete cache;
  return *this;
//EOF
_//>

//<$set_array=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::set(const double *a) {

  __CACHE_OPERATIONS__;

  if (!was_cached  ) {
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
    cache->pto = cache->to_elems.begin();
    op_count.get += cache->nelems;
    op_count.get += cache->nelems;
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

  __CACHE_OPERATIONS__;

  if (!was_cached  ) {
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
    cache->pto = cache->to_elems.begin();

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
			      const int m=0,INT_VAR_ARGS) {

  __CACHE_OPERATIONS__;

  if (!was_cached  ) {
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
    cache->line_cache_start = cache->prod_cache.begin();
    cache->nlines = nlines;
    cache->line_size = line_size;
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
      lc->starta = lc->linea.begin();
      if (!inc(findx,ndimsf)) break;
    }
    int ntot = nlines*line_size;
    op_count.get += ntot;
    op_count.put += nlines;
    __COUNT_OPER__;
  }

  double **pa,**pe,val,aux;
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
  if (!use_cache) delete cache;
  return *this;
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double FastMat2::__NAME___all(__OTHER_ARGS__) const {
  static FastMat2 retval(0);
  retval.__NAME__(*this __C__ __OTHER_ARGS_P__);
  return *retval.store;
}

//EOF
_//>
//< gen_sum_all(); //>

//<$in_place = <<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::__NAME__(__FUN_ARGS__) {

  __CACHE_OPERATIONS__;

  if (!was_cached  ) {
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
    cache->pto = cache->to_elems.begin();
  }

  double **to = cache->pto;
  double **to_end = to + cache->nelems;
  while (to<to_end) {
    __ELEM_OPERATIONS__;
  }

  if (!use_cache) delete cache;
  return *this;
}  
//EOF
_//>
//< in_place_all(); //>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
//<$prod = <<'//EOF';
FastMat2 & FastMat2::prod(const FastMat2 & A,const FastMat2 & B,const int m,INT_VAR_ARGS) {

  __CACHE_OPERATIONS__;

  if (!was_cached  ) {
    Indx ia,ib,ii;
    va_list ap;

    Indx Afdims,Bfdims,fdims;
    A.get_dims(Afdims);
    B.get_dims(Bfdims);

    // maxc:= maximum contracted index
    int nd=0,cd=0,maxc=0;
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
	     " dimes of  result.\n");
      exit(1);
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
    cache->line_cache_start = cache->prod_cache.begin();
    cache->line_size = mem_size(ndimsc);
    int line_size = cache->line_size;

    int jlc=0,inca,incb;
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
      lc->starta = lc->linea.begin();
      lc->startb = lc->lineb.begin();
      lc->inca = inca;
      lc->incb = incb;
      // lc->linear = 0; // force non-linear

      if (!inc(findx,ndimsf)) break;
    }
    // Hay que contar mejor cuantos elementos hay que traer
    // op_count.get += cache->nlines*cache->line_size;
    int ntot = cache->nlines*cache->line_size;
    op_count.put += cache->nlines*cache->line_size;
    op_count.sum += ntot;
    op_count.mult += ntot;

  }

  // Perform computations using cached addresses
  int nlines = cache->nlines;
  double **pa,**pb,**pa_end,sum,*paa,*pbb,*paa_end;
  for (int j=0; j<nlines; j++) {
    lc = cache->line_cache_start+j;
    pa = lc->starta;
    pb = lc->startb;
    if (lc->linear) {
      paa = *pa;
      pbb = *pb;
      paa_end = paa + lc->inca * cache->line_size;
      sum=0.;
      while (paa<paa_end) {
	sum += (*paa)*(*pbb);
	paa += lc->inca;
	pbb += lc->incb;
      }
    } else {
      pa_end = pa + cache->line_size;
      sum=0.;
      while (pa<pa_end) {
	sum += (**pa++)*(**pb++);
      }
    }
    *(lc->target) = sum;
  }

  if (!use_cache) delete cache;
  return *this;
}
//EOF
_//>

//< print template_subst($prod); //>//

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

  __CACHE_OPERATIONS__;

  if (!was_cached  ) {
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
    cache->pfrom = cache->from_elems.begin();
  }
  __DEFINE_TARGET_POINTER__;

  double **pfrom,**pfrom_end;
  pfrom = cache->pfrom;
  pfrom_end = pfrom + cache->nelems;
  while (pfrom < pfrom_end) {
    *to++ = **pfrom++;
  }
  if (!use_cache) delete cache;
  return *this;
}  
//EOF
_//>

//< export_vals_array();//>

//<$contraction=<<'//EOF';
FastMat2 & FastMat2::ctr(const FastMat2 & A,const int m,INT_VAR_ARGS) {

  __CACHE_OPERATIONS__;


  if (!was_cached  ) {
    Indx ia,ii;
  
    Indx Afdims,fdims;
    assert(A.defined);
    A.get_dims(Afdims);

    // maxc:= maximum contracted index
    int nd=0,cd=0,maxc=0;
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
      int nd1,nd2;
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
    cache->line_cache_start = cache->prod_cache.begin();

    int jlc=0;
    while (1) {
      lc = cache->line_cache_start + jlc++;
      lc->linea.resize(line_size);
      lc->starta = lc->linea.begin();
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
    op_count.get += ntot;
    op_count.put += cache->nlines;
    op_count.sum += ntot;
  }

  // Perform computations using cached addresses
  int nlines = cache->nlines;
  double **pa,**pa_end,sum;
  for (int j=0; j<nlines; j++) {
    lc = cache->line_cache_start+j;
    pa = lc->starta;
    pa_end = pa + cache->line_size;
    sum=0.;
    while (pa<pa_end) {
      sum += **pa++;
    }
    *(lc->target) = sum;
  }

  if (!use_cache) delete cache;
  return *this;
}
//EOF
_//>

//< print template_subst($contraction); //>//

//<$diag=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::diag(FastMat2 & A,const int m,INT_VAR_ARGS) {

  __CACHE_OPERATIONS__;

  if (!was_cached) {
    Indx ia,ii;
  
    Indx Afdims;
    assert(A.defined);
    A.get_dims(Afdims);

    int nd=0;
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
    cache->pto = cache->to_elems.begin();
    cache->pfrom = cache->from_elems.begin();
    int jj=0;
    while (1) {
      for (int k=0; k<ndims; k++) 
	aindx[k] = bindx[indxmap[k]-1];
    
      cache->to_elems[jj] = location(bindx);
      cache->from_elems[jj] = A.location(aindx);
      jj++;
      if (!inc(bindx,fdims)) break;
    }
    op_count.get += nelems;
    op_count.put += nelems;
  }

  // Perform computations using cached addresses
  double **pto, **pfrom, **pto_end;
  pto = cache->pto;
  pto_end = cache->pto + cache->nelems;
  pfrom = cache->pfrom;
  
  while  (pto<pto_end) {
    **pto++ = **pfrom++;
  }

  if (!use_cache) delete cache;
  return *this;
}
//EOF
_//>

//< print template_subst($diag); //>//

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

//<$get=<<'//EOF';
double FastMat2::get(INT_VAR_ARGS) const {

  __CACHE_OPERATIONS__;

  if (!was_cached) {
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
    op_count.get += 1;
  }

  if (!use_cache) delete cache;
  return *cache->from;
}
//EOF
_//>

//< print template_subst($get); //>//

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

//<$operator_double=<<'//EOF';
// This is somewhat redundant because you can use A.get()
FastMat2::operator double() const { 
  
  __CACHE_OPERATIONS__;

  if (!was_cached) {
    Indx fdims;
    get_dims(fdims);
    assert(fdims.size()==0);
    op_count.get += 1;
    op_count.put += 1;
  }

  if (!use_cache) delete cache;
  return *store;
}
//EOF
_//>

//< print template_subst($operator_double); //>//

//<$det=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double FastMat2::det(void) const{
  __CACHE_OPERATIONS__;

  if (!was_cached) {
    Indx dims_;
    get_dims(dims_);
    int ndims = dims_.size();
    assert (ndims==2);
    assert (dims_[0] == dims_[1]);
    int m= dims_[0];
    if (m<=3) {
      cache->nelems=m*m;
      cache->from_elems.resize(cache->nelems);
      cache->pfrom=cache->from_elems.begin();
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
  if (!use_cache) delete cache;
  return det_;
}
//EOF
_//>

//< print template_subst($det); //>//

//<$inv=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::inv(const FastMat2 & A) {

  __CACHE_OPERATIONS__;

  if (!was_cached) {
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
      cache->pfrom=cache->from_elems.begin();
      cache->to_elems.resize(cache->nelems);
      cache->pto=cache->to_elems.begin();
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
      cache->A = new Matrix(m,m);
      cache->B = new Matrix(m,m);
    }

    // This I don't know if it's correct.
    op_count.get += m*m;
    op_count.put += m*m;
    op_count.mult += m*m*m;
    op_count.sum += m*m*m;

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
    A.export_vals(*cache->A);
    *cache->B = cache->A->i();
    set(cache->B->Store());
  }

  if (!use_cache) delete cache;
  return *this;
#undef A
#undef iA
}
//EOF
_//>

//< print template_subst($inv); //>//

//<$kron=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::kron(const FastMat2 & A,const FastMat2 & B) {
  __CACHE_OPERATIONS__;
  if (!was_cached) {
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
    cache->pto = cache->to_elems.begin();
    cache->pfrom = cache->from_elems.begin();
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

    op_count.get += nelems;
    op_count.put += nelems;
    op_count.mult += nelems;
  }

  // Perform computations using cached addresses
  double **pto, **pfrom_a, **pfrom_b, **pto_end;
  pto = cache->pto;
  pto_end = cache->pto + cache->nelems;
  pfrom_a = cache->pfrom;
  pfrom_b = cache->pfrom + cache->nelems;
  
  while  (pto<pto_end) 
    (**pto++) = (**pfrom_a++) * (**pfrom_b++);

  if (!use_cache) delete cache;
  return *this;

}    
//EOF
_//>

//< print template_subst($kron); //>//


//<$eye=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::eye(const double a=1.) {
  set(0.);

  __CACHE_OPERATIONS__;

  if (!was_cached  ) {

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

    cache->pto = cache->to_elems.begin();

    op_count.put += n;
    __COUNT_OPER__;
  }

  double **pto = cache->pto;
  double **pto_end = pto + cache->nelems;

  while (pto < pto_end) {
    **pto++ = a;
  }

  if (!use_cache) delete cache;
  return *this;
}  
//EOF
_//>

//< print template_subst($eye); //>//

class detsur_cache : public FastMatSubCache {
public:
  FastMat2 g;
  int m,n;
  ~detsur_cache() {};
};

//<$detsur=<<'//EOF';
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double FastMat2::detsur() {
  __CACHE_OPERATIONS__;

  detsur_cache * dsc;
  if (!was_cached) {
    Indx fdims;
    get_dims(fdims);
    assert(fdims.size()==2);

    dsc = new detsur_cache();
    dsc->n = dim(2);
    dsc->m = dim(1);
    dsc->g.resize(2,dsc->m,dsc->m);
    cache->sc = dsc;
  }

  dsc = dynamic_cast<detsur_cache *> (cache->sc);
  if (dsc->m == 0) return 1.;
  dsc->g.prod(*this,*this,1,-1,2,-1);
  return sqrt(dsc->g.det());

}
//EOF
_//>

//< print template_subst($detsur); //>//

//<=$warn_dont_modify //>


