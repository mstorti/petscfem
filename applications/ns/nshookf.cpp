//__INSERT_LICENSE__
// $Id: nshookf.cpp,v 1.1 2003/02/07 18:51:43 mstorti Exp $
#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/dxhook.h>
#include "./nsi_tet.h"

extern GlobParam *GLOB_PARAM;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class ns_dx_hook : public dx_hook {
private:
  Time time;
public:
  Vec state() { return GLOB_PARAM->state->v(); }
  TimeData *time_data() { 
    time = GLOB_PARAM->state->t();
    return &time; 
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Hook *ns_hook_factory(const char *name) {
  Hook *hook=NULL;
  if CHECK_HOOK(ns_dx_hook);
  // else if CHECK_HOOK(dx_hook);
  return hook;
}
