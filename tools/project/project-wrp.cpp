//__INSERT_LICENSE__
// $Id: project-wrp.cpp,v 1.4 2005/10/18 02:42:25 mstorti Exp $

#include <cstdio>
#include <src/fastmat2.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <ANN/ANN.h>

#include "./project.h"
#include <libguile.h>
#include "../femref/dvector.h"
#include "../femref/guilemac.h"

scm_t_bits FemInterpTag;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
SCM_DEFINE(make_fem_interp,"make-fem-interp", 0, 0, 0,
	   (), "Creates new FEM-interp object.")
#define FUNC_NAME s_make_fem_interp
{
  FemInterp *w = new FemInterp;
  printf("created new FemInterp %p\n",w);
  SCM_RETURN_NEWSMOB(FemInterpTag,w);
}
#undef FUNC_NAME

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ "free-fem-interp"
static size_t
free_fem_interp(SCM s_w) {
  SCM_ASSERT (SCM_SMOB_PREDICATE(FemInterpTag,s_w),
              s_w, SCM_ARG1, __FUN__);
  FemInterp *w = (FemInterp *)(SCM_SMOB_DATA(s_w));
  printf("deleting FemInterp %p\n",w);
  delete w;
  return sizeof(FemInterp);
}
#undef FUNC_NAME

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
SCM_DEFINE(fem_interp_init_aux,"fem-interp-init-aux", 7, 0, 0,
	   (SCM s_w,
	    SCM s_knbr, 
	    SCM s_ndof, 
	    SCM s_ndimel,
	    SCM s_xnod,
	    SCM s_icone), 
	   "Initializes FEM-interp object from mesh data.")
#define FUNC_NAME s_fem_interp_init_aux
{
  MY_SCM_GET_ARG(w,FemInterpTag,FemInterp *,1);
  MY_SCM_GET_INT(knbr,2);
  MY_SCM_GET_INT(ndof,3);
  MY_SCM_GET_INT(ndimel,4);
  MY_SCM_GET_ARG(xnod,dvdbl_tag,dvector<double> *,5);
  MY_SCM_GET_ARG(icone,dvint_tag,dvector<int> *,6);
  // w->init(knbr,ndof,ndimel,*xnod,*icone);
  return SCM_UNSPECIFIED;
}
#undef FUNC_NAME

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
extern "C" void
init_fem_interp(void) {
  FemInterpTag = scm_make_smob_type("fem-interp",sizeof(FemInterpTag));
  scm_set_smob_free(FemInterpTag,free_fem_interp);
  scm_c_define_gsubr("make-fem-interp",0,0,0,scm_fun(make_fem_interp));
#ifndef SCM_MAGIC_SNARFER
#include "./project-wrp.x"
#endif
}
