#include <math.h>
#include <stdio.h>

#include "fem.h"
#include "fastmat2.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::get_cache_position(FastMatCachePosition & pos) {
    pos = FastMatCachePosition(cache_list,position_in_cache);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::branch(void) {

  if (!use_cache) return;
  FastMatCache *cache;
  LineCache *lc;

  if (was_cached) {
    cache = cache_list_begin[position_in_cache];
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
    position_in_cache;
#ifdef FM2_CACHE_DBG
    printf ("defining cache: ");
#endif
  }
#ifdef FM2_CACHE_DBG
  printf(" cache_list %p, cache %p, position_in_cache %d\n",
	 cache_list,cache,position_in_cache);
#endif
  
#ifdef FM2_CACHE_DBG
  printf("storing position cachelist: %p, position_in_cache: %d and branching\n",
	 cache_list,position_in_cache);
#endif
  // Store this position
  cache_list_stack.
    push_back(FastMatCachePosition(cache_list,position_in_cache));
  position_in_cache++;

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::choose(const int j) {

  FastMatCache *cache;

  if (!use_cache) return;

  cache = *(cache_list_begin+position_in_cache-1);
#ifdef FM2_CACHE_DBG
    printf("At cache %p, position_in_cache %d following branch %d\n",
	   cache,position_in_cache-1,j);
#endif

  // Pad with void pointers if needed 
  for (int jj=cache->branch.size(); jj<j+1; jj++) {
    cache->branch.push_back(NULL);
  }
  if (cache->branch[j] == NULL) {
    cache->branch[j] = new FastMatCacheList;
#ifdef FM2_CACHE_DBG
    printf("creating list %p\n",cache->branch[j]);
#endif
  }
  cache_list = cache->branch[j];

  was_cached = cache_list->list_size > 0;
  cache_list_begin = cache_list->begin();
  position_in_cache = 0;
#ifdef FM2_CACHE_DBG
    printf("was_cached: %d\n",was_cached);
#endif

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::leave(void) {

  if (!use_cache) return;
  FastMatCachePosition *last_pos;
#ifdef FM2_CACHE_DBG
  printf(" --> leaving cachelist %p",cache_list);
#endif
  last_pos = &cache_list_stack.back();
  cache_list_stack.pop_back();

  cache_list = last_pos->first;
  position_in_cache = last_pos->second+1;
  //was_cached = cache_list->size() > position_in_cache;
  was_cached = cache_list->list_size > position_in_cache;
  cache_list_begin = cache_list->begin();
#ifdef FM2_CACHE_DBG
  printf(", resuming at cache_list %p, position %d, was_cached: %d\n",
	 cache_list,position_in_cache,was_cached);
#endif
  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::resync_was_cached(void) {
  if (!use_cache) return;
  was_cached = cache_list->list_size > position_in_cache;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::jump_to(FastMatCachePosition &pos) {

  if (!use_cache) return;
#ifdef FM2_CACHE_DBG
  printf(" --> leaving cachelist %p",cache_list);
#endif

  cache_list = pos.first;
  position_in_cache = pos.second;
  //was_cached = cache_list->size() > position_in_cache;
  was_cached = cache_list->list_size > position_in_cache;
  cache_list_begin = cache_list->begin();
#ifdef FM2_CACHE_DBG
  printf(", resuming at cache_list %p, position %d, was_cached: %d\n",
	 cache_list,position_in_cache,was_cached);
#endif
  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:  
//double FastMatCacheList::operation_count(void) {
double FastMat2::operation_count(void) {
  double val=0;
  val += op_count.get;
  val += op_count.put;
  val += op_count.mult;
  val += op_count.sum;
  val += op_count.div;
  val += op_count.abs;
  val += op_count.fun;
  return val;
}

// Static members of class FastMat2 related to caches
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMatCacheList *FastMat2::cache_list = NULL;
FastMatCacheList *FastMat2::cache_list_root = NULL;
int FastMat2::use_cache = 0;
int FastMat2::was_cached = 0;
int FastMat2::was_cached_save = 0;
int FastMat2::position_in_cache=0;
FastMatCache **FastMat2::cache_list_begin=NULL;
int FastMat2::cache_list_size;
vector<FastMatCachePosition> FastMat2::cache_list_stack;
OperationCount FastMat2::op_count;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::reset_cache(void) {

  if (!use_cache) {
    op_count=OperationCount();
    return;
  } else {
    cache_list = cache_list_root;
    was_cached = cache_list->list_size > 0;
    if (!was_cached) {
      // reset operation count
      op_count=OperationCount();
    }
#ifdef FM2_CACHE_DBG
    printf("cache_list->size(), cache_list->list_size %d %d\n",
	   cache_list->size(), cache_list->list_size);
#endif
    position_in_cache = 0; 
    cache_list_begin = cache_list->begin();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::print_count_statistics(void) {
  printf("Summary of operation counts:\n"
	 "     get:  %d\n"
	 "     put:  %d\n"
	 "     sum:  %d\n"
	 "     mult: %d\n"
	 "     div:  %d\n"
	 "     abs:  %d\n"
	 "     fun:  %d\n",
	 op_count.get,
	 op_count.put,
	 op_count.sum,
	 op_count.mult,
	 op_count.div,
	 op_count.abs,
	 op_count.fun);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::activate_cache(FastMatCacheList *cache_list_=NULL) {
  use_cache=1;
  assert(cache_list_root != NULL || cache_list_ !=NULL);
  if (cache_list_!=NULL) {
    cache_list = cache_list_root = cache_list_;
    position_in_cache = 0;
    cache_list_begin = cache_list->begin();
  } else {
    was_cached=was_cached_save;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void purge_cache_list(FastMatCacheList *cache_list) {
  FastMatCache **cache_list_begin,*cache;
  FastMatCacheList cl;
  if (!cache_list) return;
#ifdef FM2_CACHE_DBG
  printf(" ---> purging cache_list %p\n",cache_list);
#endif
  int size=cache_list->size();
  cache_list_begin = cache_list->begin();
  for (int j=0; j<size; j++) {
    cache = *(cache_list_begin + j);
#ifdef FM2_CACHE_DBG
    printf(" -> purging cache %p, position %d\n",cache,j);
#endif
    for (int k=0; k< cache->branch.size(); k++) {
      purge_cache_list(cache->branch[k]);
      delete cache->branch[k];
    }
    delete cache;
  }
#ifdef FM2_CACHE_DBG
  printf(" ---> end purging cache_list %p\n",cache_list);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::void_cache(void) { 
  if (use_cache && cache_list) {
    purge_cache_list(cache_list);
    cache_list->resize(0);
    cache_list->list_size =
      cache_list_size = 0;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_int_list(const int m,va_list v,Indx *indx) {
  int k;
  for (int j=0; j<m; j++) {
    k = va_arg(v,int);
    // printf("%d -> %d\n",j,k);
    indx->push_back(k);
  }
  va_end(v);
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Perm::Perm(const int m) : vector<int>(m) {
  for (int j=0; j<m; j++) (*this)[j] = j;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Perm::print(const char *s=NULL) const {
  if (s!=NULL) printf(" %s \n");
  for (int j=0; j<this->size(); j++) {
    printf("%d -> %d, ",j+1,(*this)[j]+1);
  }
  printf("\n");
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Indx::print(const char *s=NULL) const {
  if (s!=NULL) printf("%s ",s);
  int ndims = size();
  for (int jd=0; jd<ndims; jd++) 
    printf(" %d",(*this)[jd]);
  printf("\n");
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void IndexFilter::print(void) const  {
  printf("Max val %d\n",dim);
  int size_= size();
  if (size_>0) {
    int jf;
    for (jf=0; jf<size_; jf++)
      printf(" %d",(*this)[jf]+FASTMAT_BASE);
    printf("\n");
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::FastMat2(void) {
  defined=0;
  store=NULL;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::FastMat2(const Indx & dims_) {
  create_from_indx(dims_);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::~FastMat2(void) {
  assert(defined == (store != NULL));
  if (defined) delete[] store;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::reshape(const int ndims, INT_VAR_ARGS) {

  if (was_cached) return *this;
  assert(ndims>0);
  Indx dims_;
#ifdef USE_VAR_ARGS
  va_list ap;
  va_start(ap,ndims);
  read_int_list(ndims,ap,&dims_);
#else
  READ_INT_ARG_LIST(dims_);
  assert(dims_.size()==ndims);
#endif
  int newstorage = 1;
  for (int j=0; j<ndims; j++) 
    newstorage *= dims_[j];
  assert(storage == newstorage);
  perm = Perm(ndims);
  set_indx = Indx(ndims,0);
  dims.resize(ndims);
  dims_p = dims.begin();
  n_dims = ndims;
  IndexFilter pp;
  for (int jd=0; jd<ndims; jd++) {
    int m = dims_[jd];
    assert(m>0);
    dims[jd] = pp;
    dims_p[jd].dim = m;
    perm[jd] = jd;
  }
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::resize(const int ndims, INT_VAR_ARGS) {

  // This can't be cached
  if (was_cached) {
    printf("fastmat2: can't call resize() while in cached mode\n");
    exit(0);
  }

  if (defined) {
    delete[] store; // Esta tiene problemas
    defined=0;
  }
  Indx dims_;
  // assert(ndims>0);
#ifdef USE_VAR_ARGS  
  va_list ap;
  va_start(ap,ndims);
  read_int_list(ndims,ap,&dims_);
#else
  READ_INT_ARG_LIST(dims_);
  assert(dims_.size()==ndims);
#endif

  create_from_indx(dims_);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::create_from_indx(const Indx & dims_) {
    
  int ndims = dims_.size();
  perm = Perm(ndims);
  set_indx = Indx(ndims,0);
  dims.resize(ndims);
  IndexFilter pp;
  for (int jd=0; jd<ndims; jd++) {
    int m = dims_[jd];
    assert(m>=0);
    dims[jd] = pp;
    dims[jd].dim = m;
    perm[jd] = jd;
  }
  n_dims = ndims;
  dims_p = dims.begin();
  define_matrix();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::FastMat2(const int m,INT_VAR_ARGS) {
  //  assert(m>0);
  Indx indx;
#ifdef USE_VAR_ARGS
  va_list v;
  va_start(v,m);
  read_int_list(m,v,&indx);
#else
  READ_INT_ARG_LIST(indx);
  // assert(m==indx.size());
#endif
  create_from_indx(indx);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::is(const int index,const int start=0,const int finish=0,
	       const int step=1) {
  if (was_cached) return *this;
  if (start==0) {
    // If start ==0 resets the filter
    dims[index-1].reset();
    return *this;
  }
  int dim_= dims[index-1].dim;
  int f_ = (finish>0? finish : finish<0 ? dim_ : start);
  for (int jj=start; (step>0 ? jj<=f_ : jj>=f_) ; jj+=step) 
    // Store 0 based indices 
    if (jj<=dim_)
      dims[index-1].push_back(jj-1);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::ir(const int indx,const int j=0) {
  if (was_cached) return *this;
  assert(defined);
  set_indx[indx-1] = j;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::d(const int j1,const int j2) {
  if (was_cached) return *this;
  assert(defined);
  assert(1<=j1 && j1<=n_dims);
  assert(1<=j2 && j2<=n_dims);
  assert(j1!=j2);
  assert (dims[j1-1].dim == dims[j2-1].dim);
  assert(set_indx[j1-1]==0);
  set_indx[j2-1]=-j1;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::define_matrix(void) {
  DimArray::iterator j;
  storage=1;
  for (j=dims.begin(); j!=dims.end(); j++) 
    storage *= j->dim;
  if (storage>0) {
    store = new double[storage];
    defined=1;
  } else {
    store=NULL;
    defined=0;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::print(char *s=NULL) {
  if (s!=NULL) printf("--- %s ----\n",s);
  Indx fdims;
  get_dims(fdims);
  int nd = fdims.size();
  Indx indxp;
  if (nd>2) {
    Indx ndprev(nd-2,0);
    copy(&fdims[0],&fdims[nd-2],&ndprev[0]);
    indxp = Indx(nd-2,1);
    while (1) {
      indxp.print("for first indices -> ");
      print2(indxp,fdims);
      if (!inc(indxp,ndprev)) break;
    }
  } else if  (nd==2) {
    Indx indxp;
    print2(indxp,fdims);
  } else if (nd==1) {
    Indx indxp;
    print1(indxp,fdims);
  } else if (nd==0) {
    if (!defined) {
      printf("not defined!\n");
    } else {
      printf("%f\n",*store);
    }
  }
  // return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::printd(char *s=NULL) {
  Indx fdims;
  get_dims(fdims);
  fdims.print(s);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:       
void FastMat2::print2(const Indx & indxp,const Indx & fdims) {
  int nd=fdims.size();
  Indx indx = indxp;
  indx.push_back(0);
  indx.push_back(0);
  for (int j=1; j<=fdims[nd-2]; j++) {
    for (int k=1; k<=fdims[nd-1]; k++) {
      indx[nd-2]=j;
      indx[nd-1]=k;
      // indx.print("indx: ");
      printf(" %g",*location(indx));
    }
    printf("\n");
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:     
void FastMat2::print1(const Indx & indxp,const Indx & fdims) {
  Indx indx;
  indx = indxp;
  indx.push_back(0);
  int ndims = fdims.size();
  int n = fdims[ndims-1];
  for (int j=0; j<n; j++) {
    indx[ndims-1] = j+1;
    printf(" %g",*location(indx));
  }
  printf("\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::rs() {
  if (was_cached) return *this;
  int ndims = dims.size();
  for (int jd=0; jd<ndims; jd++) {
    dims[jd].reset();
  }
  set_indx = Indx(ndims,0);
  perm = Perm(ndims);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double *FastMat2::location(const Indx & indx) const {
  int ndims = n_dims;
  // int ndims = dims.size();
  // Indx aindx(ndims,0);
  int kd=0, pos = 0,m;
  for (int jd=0; jd<ndims; jd++) {
    int l = set_indx[jd];
    if (l==0) {
      l = dims_p[jd].abs_indx(indx[perm[kd]]);
      kd++;
    } else if (l<0) {
      l=-l;
      // Compute free index
      int kdd=0;
      for (int jdd=0; jdd<l; jdd++) {
	if (set_indx[jdd]==0) kdd++;
      }
      kdd--;
      l = dims_p[l-1].abs_indx(indx[perm[kdd]]);
    } else {
      l--;
    }
    m = dims_p[jd].dim;
    assert(l>=0);
    assert(l<m);
    pos = pos*(dims_p[jd].dim) + l;
  }
  return &store[pos];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int IndexFilter::abs_indx(const int j) const {
  if (size()==0) {
    return j-FASTMAT_BASE;
  } else { 
    return (*this)[j-FASTMAT_BASE];
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::exc(const int i1,const int i2) {
  if (was_cached) return *this;
  int x = perm[i1-1];
  perm[i1-1] = perm[i2-1];
  perm[i2-1] = x;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::get_dims(Indx & fdims) const {
  int ndims = dims.size();
  fdims = Indx();
  int kd=0;
  for (int jd=0; jd<ndims; jd++) {
    if (set_indx[jd]==0) {
      kd++;
      int dim_ = dims[perm[jd]].size();
      fdims.push_back(dim_>0 ? dim_ : dims[perm[jd]].dim);
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:   
int inc(Indx & indx,const Indx &dims) {
  int ndims = dims.size();
  for (int jd=ndims-1; jd>=0; jd--) {
    indx[jd]++;
    if (indx[jd]<=dims[jd]) return 1;
    indx[jd]=1;
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FastMat2::dim(const int jd) const {
  Indx fdims; 
  get_dims(fdims); 
  assert(jd<=fdims.size());
  return fdims[jd-1];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double * FastMat2::storage_begin() {
  return store;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMatCache::FastMatCache() {

  pto=NULL;
  pfrom=NULL;
  to=NULL;
  from=NULL;
  line_cache_start=NULL;
  A=NULL;
  B=NULL;

  nelems=0;
  nlines=0;
  line_size=0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMatCache::~FastMatCache() {
  if (A) delete A;
  if (B) delete B;
}



//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
double FastMat2::sum_square_all() const {
  static FastMat2 retval(0);
  retval.sum_square(*this);
  return *retval.store;
}
#endif
