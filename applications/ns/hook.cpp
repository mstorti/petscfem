//__INSERT_LICENSE__
//$Id: hook.cpp,v 1.4 2002/09/24 02:50:48 mstorti Exp $

#ifdef USE_DLEF
#include <dlfcn.h>
#endif

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>

#include "./nsi_tet.h"

extern int MY_RANK,SIZE;

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

void dl_generic_hook::init(Mesh &mesh,Dofmap &dofmap,
			   const char *name_a) {
  name = string(name_a);
  TextHashTableFilter tf(mesh.global_options);
  tf.push(name.c_str());
  const char *filename;
  tf.get_entry("filename",filename);
  PETSCFEM_ASSERT(filename,"Couldn't find filename entry for "
		   "dl_generic_hook \"%s\"\n",name_a);  
  // Get `dlopen()' handle to the extension function
  void *handle = dlopen (filename,RTLD_LAZY);
  PETSCFEM_ASSERT(handle,"Can't dl-open file \"%s\"",
		  filename);  

  string s;
  const char *error;

#define GET_FUN(FunType,fun)					\
  s = string(name) + string("_" #fun);				\
  fun = (FunType *) dlsym(handle,s.c_str());			\
  error = dlerror();						\
  PETSCFEM_ASSERT(!error,"can't dlsym() \"%s\" in file \"%s\""	\
		  "error \"%s\"\n",s.c_str(),filename,error);  

  GET_FUN(InitFun,init_fun);
  GET_FUN(TimeStepPostFun,time_step_post_fun);
  GET_FUN(TimeStepPreFun,time_step_pre_fun);

  (*init_fun)(mesh,dofmap,name_a,fun_data);
}

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
    if (!strcmp(token,"rosi_hook")) hook = new rosi_hook;
    else PetscPrintf(PETSC_COMM_WORLD,
		     "Unknown hook \"%s\nLine: \"%s\"\n",
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

