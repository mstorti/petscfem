//__INSERT_LICENSE__
// $Id: project-wrp.cpp,v 1.1 2005/03/02 23:40:41 mstorti Exp $

#include <cstdio>
#include <src/fastmat2.h>
#include <src/dvector.h>
#include <src/dvector2.h>

#include "./project.h"
#include <libguile.h>

scm_t_bits FemInterpTag;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
SCM_DEFINE(make_fem_interp,"make-fem-interp", 0, 0, 0,
	   (), "Creates new FEM-interp object.")
#define FUNC_NAME s_make_fem_interp
{
  FemInterp *w = new FemInterp;
  SCM_RETURN_NEWSMOB(FemInterpTag,w);
}
#undef FUNC_NAME

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
SCM_DEFINE(free_fem_interp,"free-fem-interp", 1, 0, 0,
	   (SCM s_w),
	   "Frees fem interp object.")
#define FUNC_NAME s_free_fem_interp
{
  SCM_ASSERT (SCM_SMOB_PREDICATE(FemInterpTag,s_w),
              s_w, SCM_ARG1, __FUN__);
  FemInterpTag *w = (FemInterpTag *)(SCM_SMOB_DATA(s_w));
  delete w;
  return sizeof(FemInterp);
}
#undef FUNC_NAME

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
extern "C" void
init_fem_interp(void) {
  FemInterpTag = scm_make_smob_type("fem-interp",sizeof(FemInterpTag));
  scm_set_smob_free(FemInterpTag,free_fem_interp);
  scm_c_define_gsubr("make-fem-interp",0,0,0,scm_fun(make_fem_interp));
#ifndef SCM_MAGIC_SNARFER
#include "./gsguile.x"
#endif
}
