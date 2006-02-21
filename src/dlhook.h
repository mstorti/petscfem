// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: dlhook.h,v 1.4 2006/02/21 11:00:33 mstorti Exp $

#ifndef DLHOOK_H
#define DLHOOK_H

#ifdef USE_DLEF
#include <dlfcn.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class dl_generic_hook : public Hook {
public:
  dl_generic_hook() : 
    options(NULL), handle(NULL), 
    fun_data(NULL), init_fun(NULL),
    time_step_post_fun(NULL),
    time_step_pre_fun(NULL),
    stage_fun(NULL), close_fun(NULL) {}
  ~dl_generic_hook() { delete options; }
  typedef void InitFun(Mesh &mesh,Dofmap &dofmap,
		       const char *name,TextHashTableFilter *options,
		       void *&fun_data);
  typedef void TimeStepPostFun(double time,int step,
			       const vector<double> &gather_values,
			       void *fun_data);
  typedef void TimeStepPreFun(double time,int step,
			      void *fun_data);
  typedef void StageFun(const char *jobinfo,int stage, 
			double time,void *data,void *fun_data);
  typedef void CloseFun(void *fun_data);
private:
  void *handle;
  void *fun_data;
  InitFun *init_fun;
  TimeStepPostFun *time_step_post_fun;
  TimeStepPreFun *time_step_pre_fun;
  StageFun *stage_fun;
  CloseFun *close_fun;
protected:
  string name;
  TextHashTableFilter *options;
public:
  void init(Mesh &mesh,Dofmap &dofmap,const char *name);
  void time_step_pre(double time,int step) {
    (*time_step_pre_fun)(time,step,fun_data);
  }
  void time_step_post(double time,int step,
		      const vector<double> &gather_values) {
    (*time_step_post_fun)(time,step,gather_values,fun_data);
  }
  void stage(const char *jobinfo,
	     int stage,double time,void *data) {
    (*stage_fun)(jobinfo,stage,time,data,fun_data);
  }
  void close() { (*close_fun)(fun_data); }
};

#define DL_GENERIC_HOOK(prefix)						\
extern "C"								\
void prefix##_init_fun(Mesh &mesh,Dofmap &dofmap,			\
		       const char *name,TextHashTableFilter *options,	\
		       void *&fun_data) {				\
  fun_data = new prefix;						\
  ((prefix *)fun_data)->init(mesh,dofmap,options,name);			\
}									\
									\
extern "C" void								\
prefix##_time_step_pre_fun(double time,int step,void *fun_data) {	\
  ((prefix *)fun_data)->time_step_pre(time,step);			\
}									\
									\
extern "C" void								\
prefix##_time_step_post_fun(double time,int step,			\
			    const vector<double> &gather_values,	\
			    void *fun_data) {				\
  ((prefix *)fun_data)							\
    ->time_step_post(time,step,gather_values);				\
}									\
									\
extern "C" void								\
prefix##_stage_fun(const char *jobinfo,int stage,			\
                   double time,void *data,void *fun_data) {		\
  ((prefix *)fun_data)->stage(jobinfo,stage,time,data);			\
}									\
									\
extern "C" void								\
prefix##_close_fun(void *fun_data) { ((prefix *)fun_data)->close(); }

#endif
#endif
