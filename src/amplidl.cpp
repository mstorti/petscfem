//__INSERT_LICENSE__
//$Id: amplidl.cpp,v 1.8 2002/02/10 23:30:55 mstorti Exp $

#ifdef USE_DLEF

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
	      "---generic (dynamically loaded)\" fixa amplitude, \n"
	      "  file: \"%s\", name: \"%s\"\n"
	      "options table: \n",ext_filename.c_str(),
	      function_name.c_str());
  thash->print();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DLGeneric::init(TextHashTable *thash) {
  char *error;
  string name,s;
  FunTable::iterator f;
  int ierr;

  //o Filename of extension function
  TGETOPTDEF_S_ND(thash,string,ext_filename,default);
  assert(ext_filename!="default");

  //o Name of extension function
  TGETOPTDEF_S_ND(thash,string,function_name,_default);
  if (function_name=="_default") function_name = "";

  // Look for the filename in the table
  FileHandleTable::iterator h = file_handle_table.find(ext_filename);
  if (h!=file_handle_table.end()) {
    // If found, the get the file handle from the table 
    fh = h->second;
    handle = fh.handle;
  } else {
    // Get `dlopen()' handle to the extension function
    handle = dlopen (ext_filename.c_str(),RTLD_LAZY);
    assert(handle);
    // Build FileHandle
    fh.handle = handle;
    fh.fun_table = new FunTable;
    // Insert FileHandle in table
    file_handle_table[ext_filename] = fh;
  }

  // function_name:= name:= is the name as entered. If not empty, then
  // we append an underscore to get `name' For instance, if user does
  // not enter `function_name' entry then we look for functions
  // `init_fun', `clear_fun' and `eval_fun'. If user enters
  // `function_name my_fun' then we look for `my_fun_init_fun',
  // `my_fun_clear_fun' and `my_fun_eval_fun'
  name = function_name;
  if (function_name!="") name += "_";

  // Look for function name if already loaded
  // in table for this file
  f = fh.fun_table->find(function_name);
  if (f != fh.fun_table->end()) {
    // If it already exists, then get function
    // pointer from handle
    fuh = f->second;
    eval_fun = fuh.eval_fun;
    init_fun = fuh.init_fun;
    clear_fun = fuh.clear_fun;
  } else {
    // If it was not loaded then get pointers with `dlsym()' and
    // insert in table
    s = name + string("eval_fun");
    eval_fun = (EvalFun *) dlsym(handle,s.c_str());
    // if `eval_fun' is not defined then give error
    if ((error = dlerror()) != NULL)  {
      fputs(error, stderr);
      printf("\n");
      exit(1);
    }

    // Look for `init_fun'
    s = name + string("init_fun");
    init_fun = (InitFun *) dlsym(handle,s.c_str());

    // Look for `clear_fun'
    s = name + string("clear_fun");
    clear_fun = (ClearFun *) dlsym(handle,s.c_str());

    // Set pointers in handle and insert handle in table
    fuh.eval_fun  = eval_fun ;
    fuh.init_fun  = init_fun ;
    fuh.clear_fun = clear_fun;

    (*fh.fun_table)[function_name] = fuh;
  }

  // Call init function
  if (init_fun) (*init_fun)(thash,fun_data);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
DLGeneric::~DLGeneric() {
  // Delete option table
  delete thash;
  // call clear function if defined by user
  if (clear_fun) (*clear_fun)(fun_data);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double DLGeneric::eval(const TimeData *time_data) {
  double t = double(* (const Time *) time_data);
  // call eval function defined by user
  return (*eval_fun)(t,fun_data);
}

#endif
