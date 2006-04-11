// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: rhhook.h,v 1.1 2006/04/11 12:36:25 mstorti Exp $

#ifndef RHHOOK_H
#define RHHOOK_H

#include <src/dvector.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This hooks reads H fields from a file or
    at init time or each time step */ 
class read_hfields_hook : public Hook {
public:
  read_hfields_hook();
  ~read_hfields_hook();
  void init(Mesh &mesh,Dofmap &dofmap,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
  virtual Vec state()=0;
  virtual TimeData *time_data()=0;
};

#endif
