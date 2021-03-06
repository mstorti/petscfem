//__INSERT_LICENSE__
//$Id: amplidl.cpp,v 1.12.104.1 2007/02/19 20:23:56 mstorti Exp $

#ifdef USE_DLEF

#include <math.h>

#include "fem.h"
#include "getprop.h"
#include "dofmap.h"
#include "ampli.h"
#include "utils.h"
#include "util2.h"
#include <dlfcn.h>

using namespace std;

DLGeneric::FileHandleTable DLGeneric::file_handle_table;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DLGeneric::print() const {
  PetscPrintf(PETSCFEM_COMM_WORLD,
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
    error = dlerror();
    PETSCFEM_ASSERT(!error && handle,
		    "Can't dlopen() file \"%s\".\n"
		    "  dlerror(): \"%s\"\n",
		    ext_filename.c_str(),error);  
    
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

#define DL_FUN_ASSERT								\
    error = dlerror();								\
    PETSCFEM_ASSERT(!error,							\
		    "Can't load dynamic function \"%s\" in file \"%s\".\n"	\
		    "  dlerror(): \"%s\"\n",					\
		    s.c_str(),ext_filename.c_str(),error)

    // If it was not loaded then get pointers with `dlsym()' and
    // insert in table
    s = name + string("eval_fun");
    eval_fun = (EvalFun *) dlsym(handle,s.c_str());
    // if `eval_fun' is not defined then give error
    DL_FUN_ASSERT;

    // Look for `init_fun'
    s = name + string("init_fun");
    init_fun = (InitFun *) dlsym(handle,s.c_str());
    DL_FUN_ASSERT;

    // Look for `clear_fun'
    s = name + string("clear_fun");
    clear_fun = (ClearFun *) dlsym(handle,s.c_str());
    DL_FUN_ASSERT;

    // Set pointers in handle and insert handle in table
    fuh.eval_fun  = eval_fun ;
    fuh.init_fun  = init_fun ;
    fuh.clear_fun = clear_fun;

    (*fh.fun_table)[function_name] = fuh;
  }

  // Call init function
  fun_data = this;
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
double DLGeneric::eval(const TimeData *time_data,int node_a, int field_a) {
  double t = double(* (const Time *) time_data);
  node_m = node_a;
  field_m = field_a;
  // call eval function defined by user
  return (*eval_fun)(t,fun_data);
}

#endif
