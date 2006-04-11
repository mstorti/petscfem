//__INSERT_LICENSE__
// $Id: nshookf.cpp,v 1.2 2006/04/11 13:23:00 mstorti Exp $
#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/dxhook.h>
#include <src/rhhook.h>
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
class ns_read_hfields_hook : public read_hfields_hook {
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
  else if CHECK_HOOK(ns_read_hfields_hook);
  return hook;
}
