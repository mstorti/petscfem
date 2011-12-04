//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(__MATS__
               const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
    __CTX_CHECK__

    Indx indx;
    indx.push_back(m);
    // READ_INT_ARG_LIST(indx);
    READ_ARG_LIST(arg,indx,INT_ARG_LIST_DEFAULT_VAL,EXIT2);
    ctx->check(indx);
  }
#endif

  FastMatCache *cache = ctx->step();
  multiprod_subcache_t *mpsc=NULL;
  if (!ctx->was_cached) {
    
    assert(!cache->sc);
    vector<const FastMat2 *> mat_list;
    __PACK__
    
    Indx indx;
    indx.push_back(m);
    READ_INT_ARG_LIST(indx);
    mpsc = new multiprod_subcache_t(cache,ctx,*this,mat_list,indx);
    cache->sc = mpsc;
  }

  mpsc = dynamic_cast<multiprod_subcache_t *> (cache->sc);
  assert(mpsc);
  mpsc->make_prod();

  if (!ctx->use_cache) delete cache;

  return *this;
}
