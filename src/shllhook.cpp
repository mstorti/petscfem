//__INSERT_LICENSE__
//$Id: shllhook.cpp,v 1.5 2003/11/24 03:40:58 mstorti Exp $

#include <string>
#include <cstdlib>

#include <src/fem.h>
#include <src/shllhook.h>

extern int MY_RANK, SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void shell_hook::init(Mesh &mesh,Dofmap &dofmap,const char *name) {
  if (!MY_RANK) {
    AutoString commi;
    commi.sprintf("make --no-print-directory "
		  "petscfem_step=%%2$d petscfem_time=%%3$f %s_%%1$s",name);
    string comm(commi.str());
    get_string(mesh.global_options,name,comm,1);
    // PETSCFEM_ASSERT(comm,"shell_hook: Couldn't find command for hook %s\n",name);
    command_pattern.cat(comm.c_str());
    command.sprintf(command_pattern.str(),"init",-1,0.);
    system(command.str());
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void shell_hook::time_step_pre(double time,int step) {
  if (!MY_RANK) {
    command.sprintf(command_pattern.str(),"pre",step,time);
    system(command.str());
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void shell_hook::time_step_post(double time,int step,
				const vector<double> &gather_values) {
  if (!MY_RANK) {
    command.sprintf(command_pattern.str(),"post",step,time);
    system(command.str());
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void shell_hook::close() {
  if (!MY_RANK) {
    command.sprintf(command_pattern.str(),"close",-2,0.);
    system(command.str());
    command.clear();
    command_pattern.clear();
  }
}
