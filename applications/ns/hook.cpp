//__INSERT_LICENSE__
//$Id: hook.cpp,v 1.5 2002/09/24 17:59:56 mstorti Exp $

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>

#include "./nsi_tet.h"
#include "./dlhook.h"

extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void HookList::init(Mesh &mesh,Dofmap &dofmap) {
  Hook *hook;
  const char *line;
  mesh.global_options->get_entry("hook_list",line);
  if (!line) return;
  char *lcpy = local_copy(line);
  int n=0; 
  while (1) { 
    char *token = strtok((n++ == 0 ? lcpy : NULL),"[ \t\n]");
    if (!token) break;

#define CHECK_HOOK(hook_name) 				\
  (!strcmp(token,#hook_name)) hook = new hook_name 

    if CHECK_HOOK(dl_generic_hook);
    // else if CHECK_HOOK(dl_generic_hook);
    else PETSCFEM_ERROR("Unknown hook \"%s\nLine: \"%s\"\n",
			token,line);
    PETSCFEM_ASSERT(hook,"Couldn't create hook \"%s\"\n",token);

    token = strtok((n++ == 0 ? lcpy : NULL),"[ \t\n]");
    hook->init(mesh,dofmap,token);
    push_back(hook);
  }
  delete[] lcpy;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void HookList::time_step_pre(double time,int step) {
  HookList::iterator q;
  for (q=begin(); q!=end(); q++) 
    (*q)->time_step_pre(time,step);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void HookList::time_step_post(double time,int step,
			      const vector<double> &gather_values) {
  HookList::iterator q;
  for (q=begin(); q!=end(); q++) 
    (*q)->time_step_post(time,step,gather_values);
}

