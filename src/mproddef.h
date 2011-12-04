// This file was automatically created by mprod.pl
// It is supposed not to be modified by hand. 
// It is not a standard header file. You should
// not include it in a program, this is only
// included once in file fastmat2c.cpp


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,

               const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
        ctx->check(&A0);
    ctx->check(&A1);


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
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,

               const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
        ctx->check(&A0);
    ctx->check(&A1);
    ctx->check(&A2);


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
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,
                const FastMat2 & A3,

               const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
        ctx->check(&A0);
    ctx->check(&A1);
    ctx->check(&A2);
    ctx->check(&A3);


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
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,
                const FastMat2 & A3,
                const FastMat2 & A4,

               const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
        ctx->check(&A0);
    ctx->check(&A1);
    ctx->check(&A2);
    ctx->check(&A3);
    ctx->check(&A4);


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
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);
  mat_list.push_back(&A4);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
FastMat2 & 
FastMat2::prod(                const FastMat2 & A0,
                const FastMat2 & A1,
                const FastMat2 & A2,
                const FastMat2 & A3,
                const FastMat2 & A4,
                const FastMat2 & A5,

               const int m,INT_VAR_ARGS_ND) {

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
        ctx->check(&A0);
    ctx->check(&A1);
    ctx->check(&A2);
    ctx->check(&A3);
    ctx->check(&A4);
    ctx->check(&A5);


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
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);
  mat_list.push_back(&A4);
  mat_list.push_back(&A5);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
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

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
        ctx->check(&A0);
    ctx->check(&A1);
    ctx->check(&A2);
    ctx->check(&A3);
    ctx->check(&A4);
    ctx->check(&A5);
    ctx->check(&A6);


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
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);
  mat_list.push_back(&A4);
  mat_list.push_back(&A5);
  mat_list.push_back(&A6);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
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

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
        ctx->check(&A0);
    ctx->check(&A1);
    ctx->check(&A2);
    ctx->check(&A3);
    ctx->check(&A4);
    ctx->check(&A5);
    ctx->check(&A6);
    ctx->check(&A7);


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
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);
  mat_list.push_back(&A4);
  mat_list.push_back(&A5);
  mat_list.push_back(&A6);
  mat_list.push_back(&A7);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
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

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
        ctx->check(&A0);
    ctx->check(&A1);
    ctx->check(&A2);
    ctx->check(&A3);
    ctx->check(&A4);
    ctx->check(&A5);
    ctx->check(&A6);
    ctx->check(&A7);
    ctx->check(&A8);


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
      mat_list.push_back(&A0);
  mat_list.push_back(&A1);
  mat_list.push_back(&A2);
  mat_list.push_back(&A3);
  mat_list.push_back(&A4);
  mat_list.push_back(&A5);
  mat_list.push_back(&A6);
  mat_list.push_back(&A7);
  mat_list.push_back(&A8);

    
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


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
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

#ifndef NDEBUG
  if (ctx->do_check_labels) {
    ctx->check_clear();
    ctx->check("prod_mat_wrapper",this);
        ctx->check(&A0);
    ctx->check(&A1);
    ctx->check(&A2);
    ctx->check(&A3);
    ctx->check(&A4);
    ctx->check(&A5);
    ctx->check(&A6);
    ctx->check(&A7);
    ctx->check(&A8);
    ctx->check(&A9);


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
