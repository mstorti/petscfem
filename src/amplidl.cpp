//__INSERT_LICENSE__
//$Id: amplidl.cpp,v 1.4 2002/02/10 12:47:00 mstorti Exp $

#include <math.h>

#include "fem.h"
#include "getprop.h"
#include "dofmap.h"
#include "ampli.h"
#include "utils.h"
#include "util2.h"
#include <dlfcn.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DLGeneric::print() const {
  PetscPrintf(PETSC_COMM_WORLD,
	      "\"generic DL\" fixa amplitude, options table: \n");
  thash->print();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DLGeneric::init(TextHashTable *thash) {
  char *error;
  handle = dlopen ("/home/mstorti/PETSC/petscfem/test/aquifer/fun.efn",RTLD_LAZY);
  assert(handle);
  fun = (EvalFun *) dlsym(handle,"eval_fun");
  if ((error = dlerror()) != NULL)  {
    fputs(error, stderr);
    printf("\n");
    exit(1);
  }

  init_fun = (InitFun *) dlsym(handle,"init_fun");
#if 0
  if ((error = dlerror()) != NULL)  {
    fputs(error, stderr);
    printf("\n");
    exit(1);
  }
#endif

  clear_fun = (ClearFun *) dlsym(handle,"clear_fun");
#if 0
  if ((error = dlerror()) != NULL)  {
    fputs(error, stderr);
    printf("\n");
    exit(1);
  }
#endif

  if (init_fun) (*init_fun)(thash,fun_data);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
DLGeneric::~DLGeneric() {
  delete thash;
  if (clear_fun) (*clear_fun)(fun_data);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double DLGeneric::eval(const TimeData *time_data) {
  double t = double(* (const Time *) time_data);
  return (*fun)(t,fun_data);
}
