// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: dxhook.h,v 1.2 2003/02/04 13:32:01 mstorti Exp $

#ifndef DXHOOK_H
#define DXHOOK_H

#ifdef USE_DLEF
#include <dlfcn.h>
#endif

#ifdef USE_SSL
#include <HDR/sockets.h>
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class dx_hook : public Hook {
#ifdef USE_SSL
private:
  TextHashTableFilter *options;
  Socket *srvr_root,*srvr;
  Mesh *mesh;
public:
  dx_hook() : options(NULL), srvr_root(NULL) {}
  ~dx_hook() { delete options; }
  void init(Mesh &mesh,Dofmap &dofmap,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
#else
public:
  void init(Mesh &mesh,Dofmap &dofmap,const char *name) {
    PETSCFEM_ERROR0("Hook unavailable. Code not compiled with sockets library!!");  
    PetscFinalize();
    exit(0);
  }
#endif
};

#endif
