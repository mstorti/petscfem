// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: dlhook.h,v 1.2 2002/09/24 17:59:56 mstorti Exp $

#ifndef DLHOOK_H
#define DLHOOK_H

#ifdef USE_DLEF
#include <dlfcn.h>
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class dl_generic_hook : public Hook {
public:
  typedef void InitFun(Mesh &mesh,Dofmap &dofmap,
		      const char *name,void *&fun_data);
  typedef void TimeStepPostFun(double time,int step,
			       const vector<double> &gather_values,
			       void *fun_data);
  typedef void TimeStepPreFun(double time,int step,
			      void *fun_data);
private:
  void *handle;
  void *fun_data;
  InitFun *init_fun;
  TimeStepPostFun *time_step_post_fun;
  TimeStepPreFun *time_step_pre_fun;
protected:
  string name;
public:
  void init(Mesh &mesh,Dofmap &dofmap,const char *name);
  void time_step_pre(double time,int step) {
    (*time_step_pre_fun)(time,step,fun_data);
  }
  void time_step_post(double time,int step,
		      const vector<double> &gather_values) {
    (*time_step_post_fun)(time,step,gather_values,fun_data);
  }
};

#endif
