//__INSERT_LICENSE__
//$Id: amplidl.cpp,v 1.5 2002/02/10 15:29:12 mstorti Exp $

#include <math.h>

#include "fem.h"
#include "getprop.h"
#include "dofmap.h"
#include "ampli.h"
#include "utils.h"
#include "util2.h"
#include <dlfcn.h>

DLGeneric::FileHandleTable DLGeneric::file_handle_table;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DLGeneric::print() const {
  PetscPrintf(PETSC_COMM_WORLD,
	      "\"generic DL\" fixa amplitude, options table: \n");
  thash->print();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DLGeneric::init(TextHashTable *thash) {
  char *error;
  FileHandle fh;
  string name,s;
  FunTable::iterator f;
  FunHandle fuh;
  int ierr;

  //o Filename of extension function
  TGETOPTDEF_S(thash,string,ext_filename,default);
  assert(ext_filename!="default");

  //o Name of extension function
  TGETOPTDEF_S(thash,string,function_name,_default);
  if (function_name=="_default") function_name = "";

  FileHandleTable::iterator h = file_handle_table.find(ext_filename);
  if (h!=file_handle_table.end()) {
    fh = h->second;
  } else {
    // Get handle to the extension function
    handle = dlopen (ext_filename.c_str(),RTLD_LAZY);
    assert(handle);
    fh.handle = handle;
    fh.fun_table = new FunTable;
    file_handle_table[ext_filename] = fh;
  }

  // function_name:= clean_name:= is the name as
  // entered. If not empty, then we append
  // an underscore to get `name'
  name = function_name;
  if (function_name!="") name += "_";

  f = fh.fun_table->find(function_name);
  if (f != fh.fun_table->end()) {
    fuh = f->second;
    eval_fun = fuh.eval_fun;
    init_fun = fuh.init_fun;
    clear_fun = fuh.clear_fun;
  } else {
    s = name + string("eval_fun");
    eval_fun = (EvalFun *) dlsym(handle,s.c_str());
    if ((error = dlerror()) != NULL)  {
      fputs(error, stderr);
      printf("\n");
      exit(1);
    }

    s = name + string("init_fun");
    init_fun = (InitFun *) dlsym(handle,s.c_str());

    s = name + string("clear_fun");
    clear_fun = (ClearFun *) dlsym(handle,s.c_str());
    fuh.eval_fun  = eval_fun ;
    fuh.init_fun  = init_fun ;
    fuh.clear_fun = clear_fun;

    (*fh.fun_table)[function_name] = fuh;
  }

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
  return (*eval_fun)(t,fun_data);
}
