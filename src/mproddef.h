

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,

               const int m,INT_VAR_ARGS_ND) {

  FastMatCache *cache = ctx->step();
  mprodwrp_subcache_t *mpwrpsc=NULL;
  if (!ctx->was_cached) {
    mpwrpsc = new mprodwrp_subcache_t;
    assert(!cache->sc);
    cache->sc = mpwrpsc;
    vector<const FastMat2 *> &mat_list = mpwrpsc->mat_list;
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,
                const FastMat2 & A3,

               const int m,INT_VAR_ARGS_ND) {

  FastMatCache *cache = ctx->step();
  mprodwrp_subcache_t *mpwrpsc=NULL;
  if (!ctx->was_cached) {
    mpwrpsc = new mprodwrp_subcache_t;
    assert(!cache->sc);
    cache->sc = mpwrpsc;
    vector<const FastMat2 *> &mat_list = mpwrpsc->mat_list;
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,
                const FastMat2 & A3,
                const FastMat2 & A4,

               const int m,INT_VAR_ARGS_ND) {

  FastMatCache *cache = ctx->step();
  mprodwrp_subcache_t *mpwrpsc=NULL;
  if (!ctx->was_cached) {
    mpwrpsc = new mprodwrp_subcache_t;
    assert(!cache->sc);
    cache->sc = mpwrpsc;
    vector<const FastMat2 *> &mat_list = mpwrpsc->mat_list;
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);
  mat_list.push_back(&A4);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,
                const FastMat2 & A3,
                const FastMat2 & A4,
                const FastMat2 & A5,

               const int m,INT_VAR_ARGS_ND) {

  FastMatCache *cache = ctx->step();
  mprodwrp_subcache_t *mpwrpsc=NULL;
  if (!ctx->was_cached) {
    mpwrpsc = new mprodwrp_subcache_t;
    assert(!cache->sc);
    cache->sc = mpwrpsc;
    vector<const FastMat2 *> &mat_list = mpwrpsc->mat_list;
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);
  mat_list.push_back(&A4);
  mat_list.push_back(&A5);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,
                const FastMat2 & A3,
                const FastMat2 & A4,
                const FastMat2 & A5,
                const FastMat2 & A6,

               const int m,INT_VAR_ARGS_ND) {

  FastMatCache *cache = ctx->step();
  mprodwrp_subcache_t *mpwrpsc=NULL;
  if (!ctx->was_cached) {
    mpwrpsc = new mprodwrp_subcache_t;
    assert(!cache->sc);
    cache->sc = mpwrpsc;
    vector<const FastMat2 *> &mat_list = mpwrpsc->mat_list;
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);
  mat_list.push_back(&A4);
  mat_list.push_back(&A5);
  mat_list.push_back(&A6);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,
                const FastMat2 & A3,
                const FastMat2 & A4,
                const FastMat2 & A5,
                const FastMat2 & A6,
                const FastMat2 & A7,

               const int m,INT_VAR_ARGS_ND) {

  FastMatCache *cache = ctx->step();
  mprodwrp_subcache_t *mpwrpsc=NULL;
  if (!ctx->was_cached) {
    mpwrpsc = new mprodwrp_subcache_t;
    assert(!cache->sc);
    cache->sc = mpwrpsc;
    vector<const FastMat2 *> &mat_list = mpwrpsc->mat_list;
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);
  mat_list.push_back(&A4);
  mat_list.push_back(&A5);
  mat_list.push_back(&A6);
  mat_list.push_back(&A7);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,
                const FastMat2 & A3,
                const FastMat2 & A4,
                const FastMat2 & A5,
                const FastMat2 & A6,
                const FastMat2 & A7,
                const FastMat2 & A8,

               const int m,INT_VAR_ARGS_ND) {

  FastMatCache *cache = ctx->step();
  mprodwrp_subcache_t *mpwrpsc=NULL;
  if (!ctx->was_cached) {
    mpwrpsc = new mprodwrp_subcache_t;
    assert(!cache->sc);
    cache->sc = mpwrpsc;
    vector<const FastMat2 *> &mat_list = mpwrpsc->mat_list;
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);
  mat_list.push_back(&A4);
  mat_list.push_back(&A5);
  mat_list.push_back(&A6);
  mat_list.push_back(&A7);
  mat_list.push_back(&A8);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,
                const FastMat2 & A3,
                const FastMat2 & A4,
                const FastMat2 & A5,
                const FastMat2 & A6,
                const FastMat2 & A7,
                const FastMat2 & A8,
                const FastMat2 & A9,

               const int m,INT_VAR_ARGS_ND) {

  FastMatCache *cache = ctx->step();
  mprodwrp_subcache_t *mpwrpsc=NULL;
  if (!ctx->was_cached) {
    mpwrpsc = new mprodwrp_subcache_t;
    assert(!cache->sc);
    cache->sc = mpwrpsc;
    vector<const FastMat2 *> &mat_list = mpwrpsc->mat_list;
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);
  mat_list.push_back(&A4);
  mat_list.push_back(&A5);
  mat_list.push_back(&A6);
  mat_list.push_back(&A7);
  mat_list.push_back(&A8);
  mat_list.push_back(&A9);

    
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
