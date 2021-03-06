///__INSERT_LICENSE__
//$Id mstorti-v6-branch-1.0.1-5-ge86f38c Wed Sep 19 13:06:15 2007 -0300$

#include <cmath>
#include <cstdio>

using namespace std;
#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/fm2stats.h>
#include <src/fm2prod.h>

FastMat2::CacheCtx1 FastMat2::global_cache_ctx1;
FastMat2::CacheCtx2 FastMat2::global_cache_ctx2;
int FastMat2::cache_dbg=0;
int FastMat2::use_cachectx2_as_default=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::init() {
  prod2_subcache_t::init();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2::CacheCtx::~CacheCtx() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2::CacheCtx::CacheCtx() 
  : use_cache(0), was_cached(0),
    do_check_labels(0), mprod_order(mixed), 
    optimal_mprod_order_max(6) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx::check_clear() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx::check(const char *) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx::check(const FastMat2 *) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx::check(int x) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx::check(const Indx &) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx
::check(const char *,const FastMat2 *,
        const FastMat2 *,const FastMat2 *) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMatCachePosition::FastMatCachePosition() {
  first = NULL;
  second = -1;
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) 
    printf ("in constructor creating cache_position: %p\n",this);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMatCachePosition::~FastMatCachePosition() {
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) 
    printf ("in destructor deleting cache_position: %p\n",this);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMatCacheList::FastMatCacheList() { 
  list_size=0; 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMatCacheList::~FastMatCacheList() {
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) 
    printf ("in destructor deleting cache_list: %p\n",this);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::CacheCtx1::get_cache_position(FastMatCachePosition & pos) {
  pos.first = cache_list;
  pos.second = position_in_cache;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::get_cache_position(FastMatCachePosition & pos) {
  global_cache_ctx1.get_cache_position(pos);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FastMat2::CacheCtx1
::check_cache_position(FastMatCachePosition & pos) {
  if (!pos.first) {
    get_cache_position(pos);
    return 1;
  } else {
    FastMatCachePosition pos2;
    get_cache_position(pos2);
    return pos2==pos;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::CacheCtx1::deactivate_cache(void) {
  was_cached_save = was_cached;
  use_cache=0; 
  was_cached=0; 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::deactivate_cache() {
  global_cache_ctx1.deactivate_cache();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::CacheCtx1::CacheCtx1() 
  : cache_list_root(NULL), 
    cache_list(NULL), 
    position_in_cache(0),
    cache_list_begin(NULL),
    cache_list_size(0),
    was_cached_save(0),
    last_trace(""),
    do_trace(0) { use_cache=0; was_cached=0; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::branch() {
  global_cache_ctx1.branch();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::CacheCtx1::branch() {

  if (!use_cache) return;
  FastMatCache *cache;

  if (was_cached) {
    cache = cache_list_begin[position_in_cache];
#ifdef FM2_CACHE_DBG
    if (FastMat2::cache_dbg) printf ("defining cache: ");
#endif
  } else if (!use_cache) {
    cache = new FastMatCache;
  } else {
    cache = new FastMatCache;
    cache_list->push_back(cache);
    cache_list_begin = &*(cache_list->begin());
    cache_list->list_size =
      cache_list_size = cache_list->size();
#ifdef FM2_CACHE_DBG
    if (FastMat2::cache_dbg) printf ("defining cache: ");
#endif
  }
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) 
    printf(" cache_list %p, cache %p, position_in_cache %d\n",
           cache_list,cache,position_in_cache);
#endif
  
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) printf("storing position cachelist: %p, position_in_cache: %d and branching\n",
	 cache_list,position_in_cache);
#endif
  // Store this position
  cache_list_stack.push_back(FastMatCachePosition());
  FastMatCachePosition &cp = cache_list_stack.back();
  cp.first = cache_list;
  cp.second = position_in_cache;
#ifdef FM2_CACHE_DBG
    if (FastMat2::cache_dbg) 
      printf("in branch: pushing cache_pos %p, stack %p, pos %d\n",
             &cache_list_stack.back(),
             &cache_list_stack,
             cache_list_stack.size()-1);
#endif
  position_in_cache++;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::choose(const int j) {
  global_cache_ctx1.choose(j);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::CacheCtx1::choose(const int j) {

  FastMatCache *cache;

  if (!use_cache) return;

  cache = *(cache_list_begin+position_in_cache-1);
#ifdef FM2_CACHE_DBG
    if (FastMat2::cache_dbg) 
      printf("At cache %p, position_in_cache"
             " %d following branch %d\n",
             cache,position_in_cache-1,j);
#endif

  // Pad with void pointers if needed 
  for (int jj=cache->branch.size(); jj<j+1; jj++) {
    cache->branch.push_back(NULL);
  }
  if (cache->branch[j] == NULL) {
    cache->branch[j] = new FastMatCacheList;
#ifdef FM2_CACHE_DBG
    if (FastMat2::cache_dbg) 
      printf("creating list %p\n",cache->branch[j]);
#endif
  }
  cache_list = cache->branch[j];

  was_cached = cache_list->list_size > 0;
  cache_list_begin = &*(cache_list->begin());
  position_in_cache = 0;
#ifdef FM2_CACHE_DBG
    if (FastMat2::cache_dbg) 
      printf("was_cached: %d\n",was_cached);
#endif

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::leave() {
  global_cache_ctx1.leave();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::CacheCtx1::leave() {

  if (!use_cache) return;
  FastMatCachePosition *last_pos;
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) printf(" --> leaving cachelist %p",cache_list);
#endif
  last_pos = &cache_list_stack.back();

  cache_list = last_pos->first;
  position_in_cache = last_pos->second+1;
  //was_cached = cache_list->size() > position_in_cache;
  was_cached = cache_list->list_size > position_in_cache;
  cache_list_begin = &*cache_list->begin();
  cache_list_stack.pop_back();
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) printf(", resuming at cache_list %p, position %d, was_cached: %d\n",
	 cache_list,position_in_cache,was_cached);
#endif
  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::resync_was_cached(void) {
  global_cache_ctx1.resync_was_cached();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void FastMat2::CacheCtx1::resync_was_cached(void) {
  if (!use_cache) return;
  was_cached = cache_list->list_size > position_in_cache;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::jump_to(FastMatCachePosition &pos) {
  global_cache_ctx1.jump_to(pos);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::CacheCtx1::jump_to(FastMatCachePosition &pos) {

  if (!use_cache) return;
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) printf(" --> leaving cachelist %p",cache_list);
#endif

  cache_list = pos.first;
  position_in_cache = pos.second;
  //was_cached = cache_list->size() > position_in_cache;
  was_cached = cache_list->list_size > position_in_cache;
  cache_list_begin = &*cache_list->begin();
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) printf(", resuming at cache_list %p, position %d, was_cached: %d\n",
	 cache_list,position_in_cache,was_cached);
#endif
  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:  
double FastMat2::operation_count(void) {
  return global_cache_ctx1.operation_count();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:  
double FastMat2::CacheCtx1::operation_count(void) {
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::reset_cache(void) {
  global_cache_ctx1.reset_cache();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::CacheCtx1::reset_cache(void) {

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
    if (FastMat2::cache_dbg) printf("cache_list->size(), cache_list->list_size %d %d\n",
	   cache_list->size(), cache_list->list_size);
#endif
    position_in_cache = 0; 
    cache_list_begin = &*cache_list->begin();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::print_count_statistics(void) {
  global_cache_ctx1.print_count_statistics();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::CacheCtx1::print_count_statistics(void) {
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
void FastMat2
::activate_cache(FastMatCacheList *clp) {
  global_cache_ctx1.activate_cache(clp);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::CacheCtx1
::activate_cache(FastMatCacheList *cache_list_) {
  use_cache=1;
  assert(cache_list_root != NULL || cache_list_ !=NULL);
  if (cache_list_!=NULL) {
    cache_list = cache_list_root = cache_list_;
    position_in_cache = 0;
    cache_list_begin = &*cache_list->begin();
  } else {
    was_cached=was_cached_save;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::CacheCtx1::init() {
  activate_cache(&cache_list_internal);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void purge_cache_list(FastMatCacheList *cache_list) {
  FastMatCache **cache_list_begin,*cache;
  if (!cache_list) return;
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) printf(" ---> purging cache_list %p\n",cache_list);
#endif
  int size=cache_list->size();
  cache_list_begin = &*cache_list->begin();
  for (int j=0; j<size; j++) {
    cache = *(cache_list_begin + j);
#ifdef FM2_CACHE_DBG
    if (FastMat2::cache_dbg) 
      printf(" -> purging cache %p, position %d\n",cache,j);
#endif
    for (unsigned int k=0; k< cache->branch.size(); k++) {
      FastMatCacheList *cl_p =cache->branch[k];
      if (cl_p) {
#ifdef FM2_CACHE_DBG
    if (FastMat2::cache_dbg) 
      printf("deleting cache_list %p\n",cl_p);
#endif
        purge_cache_list(cl_p);
        delete cl_p;
      }
    }
    delete cache;
  }
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) printf(" ---> end purging cache_list %p\n",cache_list);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::void_cache(void) { 
  global_cache_ctx1.void_cache();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::CacheCtx1::void_cache(void) { 
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
void Perm::print(const char *s) const {
  if (s!=NULL) printf(" %s \n",s);
  for (unsigned int j=0; j<this->size(); j++) {
    printf("%d -> %d, ",j+1,(*this)[j]+1);
  }
  if (s!=NULL) printf("\n");
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::set_default_ctx(void) {
  if (use_cachectx2_as_default) 
    ctx = &global_cache_ctx2;
  else 
    ctx = &global_cache_ctx1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::FastMat2(void) 
  : store(NULL), storage_is_external(0), defined(0)  { 
  set_default_ctx();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::FastMat2(const Indx & dims_) {
  create_from_indx(dims_);
  set_default_ctx();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::FastMat2(CacheCtx *ctx_a) 
  : ctx(ctx_a), store(NULL), storage_is_external(0), defined(0)  { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::FastMat2(CacheCtx *ctx_a,const Indx & dims_) 
  : ctx(ctx_a) {
  create_from_indx(dims_);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::set_external_storage(double *ap) {
  assert(defined == (store != NULL || storage_is_external));
  assert(defined);
  if (store && !storage_is_external) delete[] store;
  store = ap;
  storage_is_external = 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::~FastMat2(void) {
  // assert(defined == (store != NULL || storage_is_external));
  if (store && !storage_is_external) delete[] store;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::reshape(const int ndims, INT_VAR_ARGS_ND) {

  if (ctx->was_cached) return *this;
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
  dims_p = &*dims.begin();
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
FastMatSubCache::~FastMatSubCache() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::resize(const int ndims, INT_VAR_ARGS_ND) {

  // This can't be cached
  if (ctx->was_cached) {
    printf("fastmat2: can't call resize() while in cached mode\n");
    abort();
    exit(0);
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

  resize(dims_);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::resize(const Indx &indx) {

  assert(!storage_is_external);
  // This can't be cached
  if (ctx->was_cached) {
    printf("fastmat2: can't call resize() while in cached mode\n");
    abort();
    exit(0);
  }

  if (defined) {
    delete[] store;
    defined=0;
  }

  create_from_indx(indx);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::clear() {

  assert(!storage_is_external);
  // This can't be cached
  if (ctx->was_cached) {
    printf("fastmat2: can't call clear() while in cached mode\n");
    exit(0);
  }

  if (defined) {
    delete[] store;
    store = NULL;
    defined=0;
  }
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
  dims_p = &*dims.begin();
  define_matrix();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::FastMat2(const int m,
                   INT_VAR_ARGS_ND) {
  set_default_ctx();
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
FastMat2::FastMat2(CacheCtx *ctx_a,const int m,INT_VAR_ARGS_ND) 
  : ctx(ctx_a) {
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
FastMat2 & FastMat2::is(const int index,const int start,const int finish,
	       const int step) {

  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("is",this);
    ctx->check(index);
    ctx->check(start);
    ctx->check(finish);
    ctx->check(step);
  }
  FastMatCache *cache = ctx->step();

  if (!ctx->was_cached) {
    assert(defined);
    if (start==0) {
      // If start ==0 resets the filter
      dims[index-1].reset();
      return *this;
    }
    int dim_= dims[index-1].dim;
    int f_ = (finish>0? finish : finish<0 ? dim_ : start);
    for (int jj=start; (step>0 ? jj<=f_ : jj>=f_) ; jj+=step) {
      // Store 0 based indices 
      assert(jj>=1 && jj<=dim_);
      dims[index-1].push_back(jj-1);
    }
  }

  if (!ctx->use_cache) delete cache;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::ir(const int indx,const int j) {

  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("ir",this);
    ctx->check(indx);
    ctx->check(j);
  }
  FastMatCache *cache = ctx->step();

  if (!ctx->was_cached) {
    assert(defined);
    set_indx[indx-1] = j;
  }

  if (!ctx->use_cache) delete cache;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::d(const int j1,const int j2) {

  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("d",this);
    ctx->check(j1);
    ctx->check(j2);
  }
  FastMatCache *cache = ctx->step();

  if (!ctx->was_cached) {
    assert(defined);
    assert(1<=j1 && j1<=n_dims);
    assert(1<=j2 && j2<=n_dims);
    assert(j1!=j2);
    assert (dims[j1-1].dim == dims[j2-1].dim);
    assert(set_indx[j1-1]==0);
    set_indx[j2-1]=-j1;
  }

  if (!ctx->use_cache) delete cache;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int FastMat2::size() const {
  Indx dims;
  get_dims(dims);
  return comp_storage_size(dims);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FastMat2::comp_storage_size(const Indx & indx) const {
  int storage=1;
  for (int j=0; j<indx.size(); j++) 
    storage *= indx[j];
  return storage;
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
    storage_is_external = 0;
  } else {
    store=NULL;
    defined=0;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::print(const char *s) const {
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
      printf("(");
      indxp.print();
      printf(" * *)\n");
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
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::print(int rowsz,
		     const char *s) const {
  if (s!=NULL) printf("--- %s ----\n",s);
  Indx fdims;
  get_dims(fdims);
  int nd = fdims.size();
  Indx indx(nd,1);
  int counter = 0;
  assert(rowsz>0);
  while (1) {
    printf(" %g",*location(indx));
    counter++;
    if (counter==rowsz) {
      printf("\n");
      counter=0;
    }
    if (!inc(indx,fdims)) break;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::dump(FILE *stream,int rowsz) const {
  Indx fdims;
  get_dims(fdims);
  int nd = fdims.size();
  Indx indx(nd,1);
  int counter = 0;
  assert(rowsz>=0);
  while (1) {
    fprintf(stream," %g",*location(indx));
    counter++;
    if (rowsz>0 && counter==rowsz) {
      fprintf(stream,"\n");
      counter=0;
    }
    if (!inc(indx,fdims)) break;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::printd(char *s) {
  Indx fdims;
  get_dims(fdims);
  fdims.print(s);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:       
void FastMat2::print2(const Indx & indxp,const Indx & fdims) const {
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
void FastMat2::print1(const Indx & indxp,const Indx & fdims) const {
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

  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("rs",this);
  }
  FastMatCache *cache = ctx->step();

  if (!ctx->was_cached) {
    int ndims = dims.size();
    for (int jd=0; jd<ndims; jd++) {
      dims[jd].reset();
    }
    set_indx = Indx(ndims,0);
    perm = Perm(ndims);
  }

  if (!ctx->use_cache) delete cache;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double *FastMat2::location(const Indx & indx) const {
  int ndims = n_dims;
  // int ndims = dims.size();
  // Indx aindx(ndims,0);
  int kd=0, pos = 0, PFUNUSED m;
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
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("exc",this);
    ctx->check(i1);
    ctx->check(i2);
  }
  FastMatCache *cache = ctx->step();

  if (!ctx->was_cached) {
    int x = perm[i1-1];
    perm[i1-1] = perm[i2-1];
    perm[i2-1] = x;
  }

  if (!ctx->use_cache) delete cache;
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
int FastMat2::n() const {
  Indx fdims; 
  get_dims(fdims); 
  return fdims.size();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double * FastMat2::storage_begin() {
  return store;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::eps_LC() {
  // Levi-Civita density tensor

  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("eps_LC",this);
  }
  FastMatCache *cache = ctx->step();

  if (!ctx->was_cached) resize(3,3,3,3);

  set(0.);
  setel(+1.,1,2,3);
  setel(+1.,2,3,1);
  setel(+1.,3,1,2);

  setel(-1.,2,1,3);
  setel(-1.,3,2,1);
  setel(-1.,1,3,2);
 
  if (!ctx->use_cache) delete cache;
  return *this;
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
  sc=NULL;

  nelems=0;
  nlines=0;
  line_size=0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMatCache::~FastMatCache() {
  // printf("deleting FM2 %p\n",this);
  delete A;
  A = NULL;
  delete B;
  B = NULL;
  delete sc;
  sc = NULL;
#if 0
  for (unsigned k=0; k<branch.size(); k++) {
    if (branch[k]) {
      delete branch[k];
      branch[k] = NULL;
    }
  }
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMatCache *FastMat2::CacheCtx1::step() {
  FastMatCache *cache;
  if (was_cached) {
    cache = cache_list_begin[position_in_cache++];
#ifdef FM2_CACHE_DBG
    if (FastMat2::cache_dbg) printf ("reusing cache: ");
#endif
  } else if (!use_cache) {
    cache = new FastMatCache;
  } else {
    cache = new FastMatCache;
    cache_list->push_back(cache);
    cache_list_begin = &*cache_list->begin();
    cache_list->list_size =
      cache_list_size = cache_list->size();
    position_in_cache++;
    cache->trace_label = last_trace;
    if (do_trace && last_trace != "")
      printf("defining cache trace label %s\n",
             last_trace.c_str());
#ifdef FM2_CACHE_DBG
    if (FastMat2::cache_dbg) printf ("defining cache: ");
#endif
  }
#ifdef FM2_CACHE_DBG
  if (FastMat2::cache_dbg) 
    printf(" cache_list %p, cache %p, position_in_cache %d\n",
           cache_list,cache,position_in_cache-1);
#endif
  if (do_trace) {
    string &s = cache->trace_label;
    if (s != "") 
      printf("passing through trace %s\n",s.c_str());
  }
  return cache;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx1::set_trace(string s) {
  last_trace = s;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx1::set_trace(const char *label) {
  last_trace = label;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::CacheCtx1::trace(int state) {
  do_trace = state;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int FastMat2::is_nan() {
  int sz = size();
  double *a = storage_begin();
  for (int j=0; j<sz; j++) 
    if (isnan(a[j])) return 1;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void FastMat2::init123() {
  Indx dims,indx;
  get_dims(dims);
  int ndims = dims.size();
  indx = Indx(ndims,1);
  while (1) {
    double val=0.;
    for (int jj=0; jj<ndims; jj++) 
      val = val*10+indx[jj];
    setel(indx,val);
    if (!inc(indx,dims)) break;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & FastMat2::rest(const FastMat2 & A ) { return minus(A); }
