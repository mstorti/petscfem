//__INSERT_LICENSE__
//$Id: dlhook.cpp,v 1.4 2002/10/18 00:44:14 mstorti Exp $

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/texthf.h>

#include "./nsi_tet.h"
#include "./dlhook.h"

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
  PETSCFEM_ASSERT(!error && handle,"Hook %s, can't dlopen() \"%s\" in file \"%s\".\n"
		  "    Error \"%s\"\n",name_a,s.c_str(),filename.c_str(),error);  

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

  (*init_fun)(mesh,dofmap,name_a,options,fun_data);
}

