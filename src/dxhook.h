// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: dxhook.h,v 1.1 2003/02/03 15:52:56 mstorti Exp $

#ifndef DLHOOK_H
#define DLHOOK_H

#ifdef USE_DLEF
#include <dlfcn.h>
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class dx_hook : public Hook {
public:
  dx_hook() : options(NULL) {}
  ~dx_hook() { delete options; }
  void init(Mesh &mesh,Dofmap &dofmap,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
};

#endif
