//__INSERT_LICENSE__
// $Id: advhookf.cpp,v 1.1.10.1 2004/03/19 18:59:57 mstorti Exp $
#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/dxhook.h>
#include "./advective.h"

extern GlobParam *GLOB_PARAM;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class advdif_dx_hook : public dx_hook {
  Vec state() { return GLOB_PARAM->x; }
  TimeData *time_data() { return GLOB_PARAM->time; }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Hook *advdif_hook_factory(const char *name) {
  Hook *hook=NULL;
  if CHECK_HOOK(advdif_dx_hook);
  // else if CHECK_HOOK(dx_hook);
  else if CHECK_HOOK(godunov_hook);
  return hook;
}
