// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: shllhook.h,v 1.2 2003/05/04 16:28:15 mstorti Exp $

#ifndef PETSCFEM_SHLLHOOK_H
#define PETSCFEM_SHLLHOOK_H

#include <src/fem.h>
#include <src/hook.h>
#include <src/autostr.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Executes a shell command. 
class shell_hook : public Hook {
private:
  /// The command to be executed
  AutoString command_pattern, command;
public:
  shell_hook() { }
  ~shell_hook() { }
  void init(Mesh &mesh,Dofmap &dofmap,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
};

#endif
