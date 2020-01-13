//__INSERT_LICENSE__
//$Id: dlhook.cpp,v 1.3 2006/02/19 01:33:15 mstorti Exp $

#ifdef USE_DLEF
#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/dlhook.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void dl_generic_hook::init(Mesh &mesh,Dofmap &dofmap,
			   const char *name_a) {
  const char *error;
  string s;

  name = string(name_a);
  options = new TextHashTableFilter(mesh.global_options);
  options->push(name.c_str());
  string filename("");
  options->get_string("filename",filename);
  PETSCFEM_ASSERT(filename!="","Couldn't find filename entry for "
		  "dl_generic_hook \"%s\"\n",name_a);  
  // Get `dlopen()' handle to the extension function
  void *handle = dlopen(filename.c_str(),RTLD_LAZY);
  error = dlerror();
  PETSCFEM_ASSERT(!error && handle,
                  "Hook %s, can't dlopen() \"%s\" in file \"%s\".\n"
		  "    Error \"%s\"\n",name_a,s.c_str(),
                  filename.c_str(),error);  

  string prefix("");
  options->get_string("prefix",prefix);
  PETSCFEM_ASSERT(prefix!="","Couldn't find prefix entry for "
		  "dl_generic_hook \"%s\", filename \n",name_a,filename.c_str());  

#define GET_FUN(FunType,fun)						\
  s = prefix + string("_" #fun);					\
  fun = (FunType *) dlsym(handle,s.c_str());				\
  error = dlerror();							\
  PETSCFEM_ASSERT(!error,						\
		  "Hook %s, can't dlsym() \"%s\" in file \"%s\".\n"	\
		  "    Error \"%s\"\n",name_a,				\
		  s.c_str(),filename.c_str(),error);  

  GET_FUN(InitFun,init_fun);
  GET_FUN(TimeStepPostFun,time_step_post_fun);
  GET_FUN(TimeStepPreFun,time_step_pre_fun);
  GET_FUN(StageFun,stage_fun);
  GET_FUN(CloseFun,close_fun);

  (*init_fun)(mesh,dofmap,name_a,options,fun_data);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
typedef Hook* (*get_hook_fun_t)();

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void dl_generic_hook2
::time_step_pre(double time,int step) {
  hookp->time_step_pre(time,step);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void dl_generic_hook2::init(Mesh &mesh,Dofmap &dofmap,
                            const char *name_a) {
  const char *error;
  string s;

  // NAME is an identifier for the dynamic hook
  // It may or may not be equal to the prefix for the hook itself
  // because there can be several instances of the same hook type
  name = string(name_a);
  options = new TextHashTableFilter(mesh.global_options);
  options->push(name.c_str());
  // Now OPTIONS is a subset of the glbal options for this hook
  string filename("");
  // The dynamic object (EFN extension) where the code for this hook is stored
  options->get_string("filename",filename);
  PETSCFEM_ASSERT(filename!="","Couldn't find filename entry for "
		  "dl_generic_hook \"%s\"\n",name_a);  
  // Get `dlopen()' handle to the extension function
  // An error is get if the extension is not found
  void *handle = dlopen(filename.c_str(),RTLD_NOW);
  error = dlerror();
  PETSCFEM_ASSERT(!error && handle,
                  "Hook %s, can't dlopen() \"%s\" in file \"%s\".\n"
		  "    Error \"%s\"\n",name_a,s.c_str(),
                  filename.c_str(),error);
  // This is the prefix for the hook type (the C++ class for the hook)
  string prefix("");
  options->get_string("prefix",prefix);
  PETSCFEM_ASSERT(prefix!="","Couldn't find prefix entry for "
		  "dl_generic_hook \"%s\", filename \n",
                  name_a,filename.c_str());
  // The hook class must be acompanied by a function that
  // returns a pointer for the class object. We will `dlopen' it
  string sfun = prefix + "_get_hook_fun";
  get_hook_fun_t fun = (get_hook_fun_t)dlsym(handle,sfun.c_str());
  error = dlerror();							\
  PETSCFEM_ASSERT(!error,						\
		  "Hook %s, can't dlsym() \"%s\" in file \"%s\".\n"	\
		  "    Error \"%s\"\n",name_a,				\
		  sfun.c_str(),filename.c_str(),error);  
  // Get a pointer to the hook object
  hookp = fun();
  // Initialize the hook object
  hookp->init(mesh,dofmap,name_a);
}

#endif
