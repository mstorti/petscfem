//__INSERT_LICENSE__
//$Id: amplidl.cpp,v 1.1 2002/02/10 00:20:09 mstorti Exp $

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
void DLGeneric::init(TextHashTable *thash_) {
  char *error;
  handle = dlopen ("./fun.efn",RTLD_LAZY);
  assert(handle);
  fun = (EvalFun *) dlsym(handle,"eval_fun");
  if ((error = dlerror()) != NULL)  {
    fputs(error, stderr);
    printf("\n");
    exit(1);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double DLGeneric::eval(const TimeData *time_data) {
  double t = double(* (const Time *) time_data);
  return (*fun)(t);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
DLGeneric::~DLGeneric() { delete thash; }

