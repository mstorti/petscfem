// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: rhhook.h,v 1.3 2006/04/11 16:44:46 mstorti Exp $

#ifndef RHHOOK_H
#define RHHOOK_H

#include <string>
#include <src/dvector.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This hooks reads H fields from a file or
    at init time or each time step */ 
class read_hfields_hook : public Hook {
private:
  int nu, ndim, nnod;
  int read_hfields_ignore_extra_nodes;
  string read_hfields_file;
  void read_hfields();
  Mesh *mesh_p;
public:
  read_hfields_hook() : mesh_p(NULL) {}
  ~read_hfields_hook() {}
  void init(Mesh &mesh,Dofmap &dofmap,
	    const char *name);
  void time_step_pre(double time,int step);
  virtual Vec state()=0;
  virtual TimeData *time_data()=0;
};

#endif
