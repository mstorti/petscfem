// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: hook.h,v 1.1 2003/01/25 17:14:54 mstorti Exp $

#ifndef HOOK_H
#define HOOK_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Hooks are executed previous and after the time step. 
class Hook {
public:
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
  void init(Mesh &mesh,Dofmap &dofmap);
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

#endif
