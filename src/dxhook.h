// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: dxhook.h,v 1.5 2003/02/09 14:50:57 mstorti Exp $

#ifndef DXHOOK_H
#define DXHOOK_H

#define USE_THREADS

#include <pthread.h>

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
  Dofmap *dofmap;
  int step_cntr, steps, ierr;
  void re_launch_connection();
  pthread_t thread;
  enum connection_state_t {
    not_launched, not_connected, connected} connection_state_m,
    connection_state_master;
  connection_state_t connection_state();
  void set_connection_state(connection_state_t s);
public:
  dx_hook() : options(NULL), srvr_root(NULL), 
    step_cntr(0), steps(0), connection_state_m(not_launched),
    connection_state_master(not_launched) {}
  ~dx_hook() { delete options; }
  void init(Mesh &mesh,Dofmap &dofmap,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
  void *wait_connection();
  virtual Vec state()=0;
  virtual TimeData *time_data()=0;
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
