//__INSERT_LICENSE__
//$Id: amplidl.cpp,v 1.2 2002/02/10 02:33:32 mstorti Exp $

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
  InitFun *init_fun = (InitFun *) dlsym(handle,"init_fun");
  if ((error = dlerror()) != NULL)  {
    fputs(error, stderr);
    printf("\n");
    exit(1);
  }
#ifdef USE_ADVANCED_DL_EFN
  (*init_fun)(thash);
#else
  (*init_fun)();
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double DLGeneric::eval(const TimeData *time_data) {
  double t = double(* (const Time *) time_data);
  return (*fun)(t);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
DLGeneric::~DLGeneric() { delete thash; }

