//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(__MATS__
               const int m,INT_VAR_ARGS_ND) {

  FastMatCache *cache = ctx->step();
  mprodwrp_subcache_t *mpwrpsc=NULL;
  if (!ctx->was_cached) {
    mpwrpsc = new mprodwrp_subcache_t;
    assert(!cache->sc);
    cache->sc = mpwrpsc;
    vector<const FastMat2 *> &mat_list = mpwrpsc->mat_list;
    __PACK__
    
    Indx &indx = mpwrpsc->indx;
    indx.push_back(m);
    READ_INT_ARG_LIST(indx);
  }

  mpwrpsc = dynamic_cast<mprodwrp_subcache_t *> (cache->sc);
  assert(mpwrpsc);

  prod(mpwrpsc->mat_list,mpwrpsc->indx);
  if (!ctx->use_cache) delete cache;

  return *this;
}
