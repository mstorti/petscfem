//__INSERT_LICENSE__
// $Id: gsguile.cpp,v 1.1 2005/01/15 13:44:59 mstorti Exp $

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

#define TRACE(j) printf("trace %d\n",j)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUN__
#define __FUN__ "getsurf"
SCM getsurf2(SCM icone_s,SCM base_s) {

  // parse args
  SCM_ASSERT(SCM_SMOB_PREDICATE(dvint_tag,icone_s),
	     icone_s, SCM_ARG1, __FUN__);
  const dvector<int> *icone_p 
    = (const dvector<int> *)SCM_SMOB_DATA (icone_s);

  int base;
  if (base_s == SCM_UNDEFINED) base = 0;
  else {
    SCM_ASSERT(SCM_INUMP(base_s),
	       base_s, SCM_ARG2, __FUN__);
    base = SCM_INUM(base_s);
  }
  getsurf(*icone_p,base);
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
  scm_c_define_gsubr("getsurf",1,1,0,scm_fun(getsurf2));
  scm_c_define_gsubr("my-dv-print",1,0,0,scm_fun(my_dv_print));
}
