// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: dxhook.h,v 1.6 2003/02/09 22:39:57 mstorti Exp $

#ifndef DXHOOK_H
#define DXHOOK_H

#ifdef USE_PTHREADS
#include <pthread.h>
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

#ifdef USE_PTHREADS
  enum connection_state_t {
    not_launched, not_connected, connected} connection_state_m,
    connection_state_master;
  void re_launch_connection();
  pthread_t thread;
  connection_state_t connection_state();
  void set_connection_state(connection_state_t s);
#endif

public:
  dx_hook();
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
