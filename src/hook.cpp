//__INSERT_LICENSE__
//$Id: hook.cpp,v 1.12.10.1 2007/02/19 20:23:56 mstorti Exp $

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/dlhook.h>
#include <src/shllhook.h>

Hook * Hook::factory(const char *name) {
  Hook *hook=NULL;
  if (0) {} // this is tricky!!
#ifdef USE_DLEF
  else if CHECK_HOOK(dl_generic_hook);
#endif
  else if CHECK_HOOK(shell_hook);
  return hook;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void HookList::init(Mesh &mesh,Dofmap &dofmap,
		    HookFactory *hf) {
  Hook *hook;
  const char *line;
  mesh.global_options->get_entry("hook_list",line);
  if (!line) return;
  char *lcpy = local_copy(line);
  int n=0; 
  const char *token;
  char *save_ptr;
  while (1) { 
    token = strtok_r((n++ == 0 ? lcpy : NULL),"[ \t\n]",&save_ptr);
    if (!token) break;
#if 1
    hook = Hook::factory(token);
    if (!hook && hf) hook = hf(token);
    PETSCFEM_ASSERT(hook,"Couldn't create hook \"%s\"\n",token);

    token = strtok_r((n++ == 0 ? lcpy : NULL),"[ \t\n]",&save_ptr);
    PETSCFEM_ASSERT(token,"Couldn't find name for hook type \"%s\"\n"
		    "  Line: \"%s\"\n",token,line);
    hook->init(mesh,dofmap,token);
    push_back(hook);
#endif
  }
  delete[] lcpy;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
HookList::~HookList() {
  for (unsigned int j=0; j<size(); j++) {
    delete (*this)[j];
    (*this)[j] = NULL;
  }
  clear();
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void HookList::stage(const char *jobinfo,int stage, 
		     double time,void *data) {
  HookList::iterator q;
  for (q=begin(); q!=end(); q++) 
    (*q)->stage(jobinfo,stage,time,data);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void HookList::close() {
  HookList::iterator q;
  for (q=begin(); q!=end(); q++) (*q)->close();
}

