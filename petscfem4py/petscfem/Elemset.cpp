#include "Elemset.h"

#include <cstring>
#include <fem.h>
#include <readmesh.h>


//PyPF::Elemset::Elemset() : Ptr() { }
//PyPF::Elemset::Elemset(const PyPF::Elemset & es) : Ptr(es) { }

PyPF::Elemset::~Elemset()
{ 
  /* delete[] this->ptr; */
}


PyPF::Elemset::Elemset(const std::string& type, 
		       const std::string& name)
  : Ptr(0)
{
  // type
  char* _type = new char[type.size()+1];
  strcpy(_type, type.c_str());
  bless_elemset(_type, this->ptr);
  if (! (this->ptr)) {
    delete[] _type;
    throw Error("unknown elemset type");
  }
  (*this)->type = _type;
  // name
  const std::string& _name =  (name.size()) ? name : (::Elemset::anon);
  (*this)->register_name(_name, (*this)->type);

  // options
  (*this)->thash = new TextHashTable;

  (*this)->ndof  = 0; 
  (*this)->epart = NULL;
  (*this)->isfat = 0;

  // int props
  (*this)->neliprops        = 0; 
  (*this)->elemiprops       = NULL; 
  (*this)->elem_iprop_names = g_hash_table_new(&g_str_hash,&g_str_equal); 
  (*this)->neliprops_add = 0; 
  (*this)->elemiprops_add = NULL; 

  // double props
  (*this)->nelprops         = 0; 
  (*this)->elemprops        = NULL; 
  (*this)->elem_prop_names  = g_hash_table_new(&g_str_hash,&g_str_equal); 
  (*this)->nelprops_add = 0; 
  (*this)->elemprops_add  = NULL; 

}

std::string PyPF::Elemset::getType()
{
  return (*this)->type;
}

std::string PyPF::Elemset::getName()
{
  return (*this)->name();
}


void PyPF::Elemset::setUp()
{
  (*this)->initialize();
}

std::string
PyPF::Elemset::getOption(const std::string& key)
{
  const char* value = NULL;
  (*this)->thash->get_entry(key.c_str(), value);
  return value;
}

void 
PyPF::Elemset::setOption(const std::string& key,
			  const std::string& value)
{
  (*this)->thash->add_entry(key.c_str(), value.c_str());
}

void PyPF::Elemset::setConnectivity(int nelem, int nel, int icone[]) 
{
  (*this)->nelem = nelem;
  (*this)->nel   = nel; 
  (*this)->icone = new int[nelem*nel];
  memcpy((*this)->icone, icone, nelem*nel*sizeof(int));
  (*this)->elem_conne = new int[nel];
  (*this)->ndof = 0;
}

void PyPF::Elemset::getConnectivity(int* nelem, int* nel, int* icone[]) 
{
  *nelem = (*this)->nelem;
  *nel   = (*this)->nel;
  *icone = (*this)->icone;
}

int PyPF::Elemset::getSize()
{
  return (*this)->nelem;
}

int PyPF::Elemset::getNDof() 
{
  return (*this)->ndof;
}

void PyPF::Elemset::setNDof(int ndof)
{
  (*this)->ndof = ndof;
}


void
PyPF::Elemset::print()
{
  (*this)->print();
}
