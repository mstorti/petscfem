//__INSERT_LICENSE__
//$Id: shllhook.cpp,v 1.1 2003/05/04 16:20:18 mstorti Exp $

#include <string>
#include <cstdlib>

#include <src/fem.h>
#include <src/shllhook.h>

extern int MY_RANK, SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void shell_hook::init(Mesh &mesh,Dofmap &dofmap,const char *name) {
  AutoString commi;
  commi.sprintf("make petscfem_step=%%2$d petscfem_time=%%3$f %s_%%1$s",name);
  string comm(commi.str());
  get_string(mesh.global_options,name,comm,1);
  // PETSCFEM_ASSERT(comm,"shell_hook: Couldn't find command for hook %s\n",name);
  command_pattern.cat(comm.c_str());
  command.sprintf(command_pattern.str(),"init",-1,0.);
  system(command.str());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void shell_hook::time_step_pre(double time,int step) {
  command.sprintf(command_pattern.str(),"pre",step,time);
  system(command.str());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void shell_hook::time_step_post(double time,int step,
				const vector<double> &gather_values) {
  command.sprintf(command_pattern.str(),"post",step,time);
  system(command.str());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void shell_hook::close() {
  command.sprintf(command_pattern.str(),"close",-1,0.);
  system(command.str());
  command.clear();
  command_pattern.clear();
}
