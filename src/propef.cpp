//__INSERT_LICENSE__
//$Id: propef.cpp,v 1.5 2003/02/09 22:39:57 mstorti Exp $

#ifdef USE_DLEF
#include <dlfcn.h>
#endif

#include <vector>
#include <set>

#include <src/fem.h>
#include <src/utils.h>
#include <src/getprop.h>
#include <src/elemset.h>
#include <src/readmesh.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "NewElemset::get_prop"
void NewElemset::get_prop(Property &prop,const char *prop_name,int n=1) const {

  // looks int the properties-per-element table
  props_hash_entry *phe = (props_hash_entry *)
    g_hash_table_lookup(elem_prop_names,(void *)prop_name);
  if(phe!=NULL) {
  // If the name is found in the per element properties table
  // (the line props in the text_hash_table) then we store the
  // position in the table and the value will be loaded with
  // 'load_props' for each element
    prop.indx = phe->position;
    prop.length = phe->width;
    return;
  } 

  // Look as a constant in the text_hash_table.
  int ierr = get_vec_double(prop_name,prop.val,0);
  if (!ierr) {
    prop.ptr = prop.val.begin();
    prop.length = prop.val.size();
    return;
  }

  // Look in the text_hash_table as a function
  string key;
  const char *entry;
  key = string("temp_fun[") + string(prop_name) + string("]");
  thash->get_entry(key.c_str(),entry);
  if (entry) {
#ifdef USE_DLEF
    char *entry_c = new char[strlen(entry)+1];
    strcpy(entry_c,entry);
    char *token = strtok(entry_c,"[ \t\n]");
    PETSCFEM_ASSERT(token,"NewElemset::get_prop(): can't parse file name"
		     "in entry \"%s\"\n",entry);  
    // Get `dlopen()' handle to the extension function
    void *handle = dlopen (token,RTLD_LAZY);
    PETSCFEM_ASSERT(handle,"NewElemset::get_prop(): can't open file \"%s\"",
		     token);  
    assert(handle);

    token = strtok(NULL,"[ \t\n]");
    PETSCFEM_ASSERT(token,"NewElemset::get_prop(): can't parse function name"
		     "in entry \"%s\"\n",entry);
    string s;
    const char *error;

#define GET_DL_FUN(FunType,fun)					\
    s = string(token) + string("_" #fun);			\
    prop.fun = (Property::FunType *) dlsym(handle,s.c_str());	\
    error = dlerror();						\
    PETSCFEM_ASSERT(!error,"NewElemset::get_prop():"		\
		     " dlsym() error \"%s\"\n",error);  

    GET_DL_FUN(InitFun,init_fun);
    GET_DL_FUN(EvalFun,eval_fun);
    GET_DL_FUN(ClearFun,clear_fun);

    if (prop.init_fun) prop.init_fun(thash,prop.fun_data);
    delete[] entry_c;
    return;
#else
    PETSCFEM_ERROR(": error: Can't use dynamically"
    " loaded fun \"%s\"\n",key.c_str());  
#endif
  }

//    PETSCFEM_ERROR("NewElemset::get_prop(): Not found property %s\n",
//  		 prop_name);  
}
