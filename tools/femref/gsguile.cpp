//__INSERT_LICENSE__
// $Id: gsguile.cpp,v 1.3 2005/01/16 14:34:15 mstorti Exp $

#include <string>
#include <list>
#include <set>
#include <ctime>
#include <unistd.h>
#include <multimap.h>
// #include <algorithm>
#include <limits.h>
#include "./hasher.h"
#include <src/fastmat2.h>
#include <libguile.h>

using namespace std;

#define VERBOSE 1

#include "./femref.h"
#include "./gtemplates.h"
#include "./dvector.h"

typedef SCM(*scm_fun)();

scm_t_bits GetSurfCtxTag;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static SCM
make_get_surf_ctx() {
  GetSurfCtx *w = new GetSurfCtx;
  SCM_RETURN_NEWSMOB(GetSurfCtxTag, w);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ "get-surf-ctx(free)"
static size_t
free_get_surf_ctx(SCM s_w) {
  SCM_ASSERT (SCM_SMOB_PREDICATE(GetSurfCtxTag, s_w),
              s_w, SCM_ARG1, __FUN__);
  GetSurfCtx *ctx = (GetSurfCtx *)(SCM_SMOB_DATA(s_w));
  delete ctx;
  return sizeof(GetSurfCtx);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ "getsurf"
SCM getsurf2(SCM ctx_s,SCM icone_s,
	     SCM surf_con_s,SCM surf_nodes_s,
	     SCM base_s,SCM verbose_s) {

  // parse args
  SCM_ASSERT(SCM_SMOB_PREDICATE(GetSurfCtxTag,ctx_s),
	     ctx_s, SCM_ARG1, __FUN__);
  GetSurfCtx *ctx 
    = (GetSurfCtx *)SCM_SMOB_DATA (ctx_s);

  SCM_ASSERT(SCM_SMOB_PREDICATE(dvint_tag,icone_s),
	     icone_s, SCM_ARG2, __FUN__);
  const dvector<int> *icone_p 
    = (const dvector<int> *)SCM_SMOB_DATA (icone_s);

  SCM_ASSERT(SCM_SMOB_PREDICATE(dvint_tag,surf_con_s),
	     surf_con_s, SCM_ARG3, __FUN__);
  dvector<int> *surf_con_p 
    = (dvector<int> *)SCM_SMOB_DATA (surf_con_s);

  SCM_ASSERT(SCM_SMOB_PREDICATE(dvint_tag,surf_nodes_s),
	     surf_nodes_s, SCM_ARG4, __FUN__);
  dvector<int> *surf_nodes_p 
    = (dvector<int> *)SCM_SMOB_DATA (surf_nodes_s);

  int base;
  if (base_s == SCM_UNDEFINED) base = 0;
  else {
    SCM_ASSERT(SCM_INUMP(base_s),
	       base_s, SCM_ARG5, __FUN__);
    base = SCM_INUM(base_s);
  }

  int verbose;
  if (verbose_s == SCM_UNDEFINED) verbose = 0;
  else {
    SCM_ASSERT(SCM_INUMP(verbose_s),
	       verbose_s, SCM_ARG6, __FUN__);
    verbose = SCM_INUM(verbose_s);
  }

  getsurf(*ctx,*icone_p,*surf_con_p,*surf_nodes_p,
	  base,verbose);
  return SCM_UNSPECIFIED;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ "my-dv-print"
SCM my_dv_print(SCM s_w) {
  SCM_ASSERT (SCM_SMOB_PREDICATE(dvdbl_tag,s_w),
              s_w, SCM_ARG1, __FUN__);
  dvector<double> *w = (dvector<double> *)SCM_SMOB_DATA (s_w);
  for (int j=0; j<w->size(); j++) {
    printf("j: %g\n",w->ref(j));
  }
  return SCM_UNSPECIFIED;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
extern "C" void
init_femref(void) {
  GetSurfCtxTag = scm_make_smob_type("get-surf-ctx",sizeof(GetSurfCtx));
  scm_set_smob_free(GetSurfCtxTag,free_get_surf_ctx);
  scm_c_define_gsubr("make-get-surf-ctx",0,0,0,scm_fun(make_get_surf_ctx));
  scm_c_define_gsubr("getsurf",4,2,0,scm_fun(getsurf2));
  scm_c_define_gsubr("my-dv-print",1,0,0,scm_fun(my_dv_print));
}
