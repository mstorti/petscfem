// $Id: Mesh.cpp,v 1.1.2.3 2006/03/06 16:56:04 rodrigop Exp $

#include "Nodedata.h"
#include "Elemset.h"
#include "Mesh.h"

#include <fem.h>
#include <readmesh.h>

PyPF::Mesh::~Mesh() 
{ }

PyPF::Mesh::Mesh() 
{ }

PyPF::OptionTable*
PyPF::Mesh::get_opt_table(bool create)
{
  OptionTable*& options = (*this)->global_options;
  if(options == NULL && create) options = new OptionTable;
  return options;
}

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


PyPF::Nodedata
PyPF::Mesh::getNodeData()
{
  /* test */
  if ((*this)->nodedata == NULL) throw Error("null Nodedata");

  Mesh::Base* mesh = *this;

  return mesh->nodedata;
}

PyPF::Elemset
PyPF::Mesh::getElemset(int i)
{
  /* test */
  if ((*this)->elemsetlist == NULL) throw Error("null Elemset list");

  Mesh::Base* mesh = *this;

  int n = da_length(mesh->elemsetlist);
  if (n == 0) throw Error("empty Elemset list");
  if (i<0 || i>=n) throw Error("index out of range");
  ::Elemset* e = *((::Elemset **) da_ref(mesh->elemsetlist, i));
  if (e == NULL) throw Error("null item in Elemset list");
  return e;
}

PyPF::Elemset
PyPF::Mesh::getElemset(const std::string& name)
{
  /* test */
  if ((*this)->elemsetlist == NULL) throw Error("null Elemset list");

  Mesh::Base* mesh = *this;

  ::Elemset* e = mesh->find(name);
  if (e == NULL) throw Error("Elemset not found in Elemset list");
  return e;
}

PyPF::Elemset
PyPF::Mesh::addElemset(const std::string& name, 
		       const std::string& type)
{

  if ((*this)->elemsetlist == NULL) {
    (*this)->elemsetlist = da_create(sizeof(Elemset::Base*));
  }

  Mesh::Base* mesh = *this;

  // create
  Elemset::Base* elemset = NULL;
  char* _type = new char[type.size()+1];
  strcpy(_type, type.c_str());
  bless_elemset(_type, elemset);
  if (elemset == NULL) {
    delete[] _type; throw Error("unknown Elemset type");
  }
  elemset->type = _type;
  // set name
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
  elemset->thash = new OptionTable;
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
  da_append(mesh->elemsetlist, &elemset);

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
  if ((*this)->elemsetlist == NULL) return 0;
  
  Mesh::Base* mesh = *this;

  return da_length(mesh->elemsetlist);
}
  
