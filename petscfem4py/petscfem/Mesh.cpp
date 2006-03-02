// $Id: Mesh.cpp,v 1.1.2.2 2006/03/02 21:37:12 rodrigop Exp $

#include "Nodedata.h"
#include "Elemset.h"
#include "Mesh.h"

#include <fem.h>
#include <readmesh.h>

PyPF::Mesh::Mesh() { }
PyPF::Mesh::~Mesh() { }

// PyPF::Mesh::Mesh() 
//   : Ptr(new ::Mesh) 
// { 
//   (*this)->nodedata = NULL;
//   (*this)->elemsetlist = NULL; 
//   (*this)->global_options = NULL;
// }

// PyPF::Mesh::~Mesh()
// { 
// #if 0
//   if ((*this)->elemsetlist == NULL) {
//     (*this)->elemsetlist = da_create(sizeof(::Elemset *));
//   }
//   if ((*this)->global_options != NULL) {
//     delete (*this)->global_options;
//     (*this)->global_options = NULL;
//   }
//   delete this->ptr;
// #endif
// }


// std::map<std::string,std::string> 
// PyPF::Mesh::getOptions() {
//   std::map<std::string,std::string> options;
//   return options;
// }

// void 
// PyPF::Mesh::setOptions(const std::map<std::string,std::string>& options) 
// {
//   typedef std::map<std::string,std::string> mapstr;
//   mapstr::const_iterator m = options.begin();
//   while (m!=options.end()) {
//     printf("%s : %s\n", m->first.c_str(), m->second.c_str());
//   }
// }

std::string
PyPF::Mesh::getOption(const std::string& key)
{
  if ((*this)->global_options == NULL) {
    throw Error("null option table");
    (*this)->global_options = new ::TextHashTable;
    (*this)->global_options->register_name("global_options");
    (*this)->global_options->set_as_global();
  }
  const char* value = NULL;
  (*this)->global_options->get_entry(key.c_str(), value);
  if (value == NULL) throw Error("option not found");
  return value;
}

void 
PyPF::Mesh::setOption(const std::string& key,
		      const std::string& value)
{
  if ((*this)->global_options == NULL) {
    throw Error("null option table");
    (*this)->global_options = new TextHashTable;
    (*this)->global_options->register_name("global_options");
    (*this)->global_options->set_as_global();
  }
  (*this)->global_options->add_entry(key.c_str(), value.c_str());
}

PyPF::Nodedata
PyPF::Mesh::getNodeData()
{
  ::Nodedata* nodedata = (*this)->nodedata;
  if (nodedata == NULL) throw Error("null nodedata");
  return nodedata;
}

PyPF::Elemset
PyPF::Mesh::getElemset(int i)
{
  Darray* elemsetlist = (*this)->elemsetlist;
  if (elemsetlist == NULL) throw Error("null elemset list");
  int n = da_length(elemsetlist);
  if (n == 0) throw Error("empty elemset list");
  if (i<0 || i>=n) throw Error("index out of range");
  ::Elemset* e = *((::Elemset **) da_ref(elemsetlist, i));
  if (e == NULL) throw Error("null item in elemset list");
  return e;
}

PyPF::Elemset
PyPF::Mesh::getElemset(const std::string& name)
{
  Darray* elemsetlist = (*this)->elemsetlist;
  if (elemsetlist == NULL) throw Error("null elemset list");
  ::Elemset* e = (*this)->find(name);
  if (e == NULL) throw Error("elemset not found in elemset list");
  return e;
}

PyPF::Elemset
PyPF::Mesh::addElemset(const std::string& name, 
		       const std::string& type)
{
  ::Elemset* elemset = NULL;
  // type
  char* _type = new char[type.size()+1];
  strcpy(_type, type.c_str());
  bless_elemset(_type, elemset);
  if (elemset == NULL) {
    delete[] _type; throw Error("unknown elemset type");
  }
  elemset->type = _type;
  // name
  const std::string& _name =  (name.size()) ? name : (::Elemset::anon);
  elemset->register_name(_name, _type);

  elemset->nelem = 0; 
  elemset->nel   = 0; 
  elemset->ndof  = 0;
  elemset->icone      = NULL;
  elemset->elem_conne = NULL;

  elemset->epart  = NULL;
  elemset->epart2 = NULL;
  elemset->isfat  = 0;

  // options table
  elemset->thash = new TextHashTable;
  elemset->elem_iprop_names = g_hash_table_new(&g_str_hash,&g_str_equal);
  elemset->elem_prop_names  = g_hash_table_new(&g_str_hash,&g_str_equal);

  elemset->thash->register_name(elemset->name());

  // double props
  elemset->nelprops         = 0; 
  elemset->elemprops        = NULL; 
  elemset->nelprops_add     = 0; 
  elemset->elemprops_add    = NULL; 

  // int props
  elemset->neliprops        = 0; 
  elemset->elemiprops       = NULL; 
  elemset->neliprops_add    = 0; 
  elemset->elemiprops_add   = NULL; 

  // add to elemset list
  Darray* elist = (*this)->elemsetlist;
  if (elist == NULL) {
    elist = (*this)->elemsetlist = da_create(sizeof(Elemset *));
  }
  da_append(elist, &elemset);

  // return added elemset
  return elemset;
}

// void 
// PyPF::Mesh::setNodedata(PyPF::Nodedata& nodedata)
// {
//   (*this)->nodedata = nodedata;
// }

// void 
// PyPF::Mesh::addElemset(PyPF::Elemset& elemset)
// {
//   if ((*this)->elemsetlist == NULL) {
//     (*this)->elemsetlist = da_create(sizeof(::Elemset *));
//   }
//   ::Elemset* e = elemset;
//   da_append((*this)->elemsetlist, &e);
  
// }

// bool
// PyPF::Mesh::hasElemset(const std::string& name)
// {
//   if ((*this)->elemsetlist == NULL) {
//     (*this)->elemsetlist = da_create(sizeof(::Elemset *));
//   }
//   ::Elemset* e = (*this)->find(name);
//   return e?true:false;
// }

// PyPF::Elemset
// PyPF::Mesh::findElemset(const std::string& name)
// {
//   if ((*this)->elemsetlist == NULL) {
//     (*this)->elemsetlist = da_create(sizeof(::Elemset *));
//   }
//   ::Elemset* e = (*this)->find(name);
//   return e;
// }

int
PyPF::Mesh::getSize()
{
  Darray* elemsetlist = (*this)->elemsetlist;
  if (elemsetlist == NULL)  throw Error("null elemset list");
  return da_length(elemsetlist);
}
  
