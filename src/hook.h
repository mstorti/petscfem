// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: hook.h,v 1.2 2003/02/04 23:28:47 mstorti Exp $

#ifndef HOOK_H
#define HOOK_H

class Hook;
typedef Hook *HookFactory(const char *name);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Hooks are executed previous and after the time step. 
class Hook {
public:
  static Hook *factory(const char *name);
  /** Initializes the hook. 
      @param mesh (input) the mesh
      @param dofmap (input) the dofmap
      @param name (input) the name of this hook */ 
  virtual void init(Mesh &mesh,Dofmap &dofmap,const char *name) {}
  /** This is executed previously to the time step. 
      @param time (input) the time of this time step
      @param step (input) the time step number */ 
  virtual void time_step_pre(double time,int step) {}
  /** This is executed after the time step. 
      @param time (input) the time of this time step
      @param step (input) the time step number 
      @param gather_values (input) the values gathered at this time step  */ 
  virtual void time_step_post(double time,int step,
			      const vector<double> &gather_values) {}
  /** This is executed after the finalization of the program */ 
  virtual void close() {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// This is a list of hooks that are executed once at a time
class HookList : public vector<Hook *> {
public:
  /** This executes the ``init'' function for each hook. 
      @param mesh (input) the mesh
      @param dofmap (input) the dofmap */ 
  void init(Mesh &mesh,Dofmap &dofmap,HookFactory *hf = NULL);
  /** This calls the ``time_step_pre'' function for each
      hook in the list.
      @param time (input) the time of this time step
      @param step (input) the time step number */ 
  void time_step_pre(double time,int step);
  /** This calls the ``time_step_pre'' function for each hook 
      in the list.
      @param time (input) the time of this time step
      @param step (input) the time step number 
      @param gather_values (input) the values gathered at this time step  */ 
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  /** This is executed after the finalization of the program */ 
  virtual void close();
  /// Dtor.
  ~HookList();
};

#define CHECK_HOOK(hook_name)				\
  (!strcmp(name,#hook_name)) hook = new hook_name 

#endif
