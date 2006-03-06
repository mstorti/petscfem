// $Id: Elemset.cpp,v 1.1.2.4 2006/03/06 16:56:04 rodrigop Exp $

#include "Elemset.h"

#include <fem.h>
#include <readmesh.h>


//PyPF::Elemset::Elemset() : Ptr() { }
//PyPF::Elemset::Elemset(const PyPF::Elemset & es) : Ptr(es) { }

PyPF::Elemset::~Elemset()
{ }

PyPF::Elemset::Elemset()
{ }

PyPF::OptionTable*
PyPF::Elemset::get_opt_table(bool create)
{
  OptionTable*& options = (*this)->thash;
  if(options == NULL && create) options = new OptionTable;
  return options;
}


PyPF::Elemset::Elemset(const std::string& type,
		       const std::string& name)
  : Ptr(NULL)
{
  // create
  Elemset::Base*& elemset = *this;
  char* _type = new char[type.size()+1];
  strcpy(_type, type.c_str());
  bless_elemset(_type, elemset);
  if (elemset == NULL) {
    delete[] _type; throw Error("unknown elemset type");
  }
  elemset->type = _type;
  const std::string& _name =  (name.size()) ? name : (Elemset::Base::anon);
  elemset->register_name(_name, elemset->type);

  // options
  elemset->thash = new OptionTable;

  elemset->ndof  = 0; 
  elemset->epart = NULL;
  elemset->isfat = 0;

  // int props
  elemset->neliprops        = 0; 
  elemset->elemiprops       = NULL; 
  elemset->elem_iprop_names = g_hash_table_new(&g_str_hash,&g_str_equal); 
  elemset->neliprops_add = 0; 
  elemset->elemiprops_add = NULL; 

  // double props
  elemset->nelprops         = 0; 
  elemset->elemprops        = NULL; 
  elemset->elem_prop_names  = g_hash_table_new(&g_str_hash,&g_str_equal); 
  elemset->nelprops_add = 0; 
  elemset->elemprops_add  = NULL; 

}

std::string 
PyPF::Elemset::getType()
{
  return (*this)->type;
}

std::string 
PyPF::Elemset::getName()
{
  return (*this)->name();
}


void 
PyPF::Elemset::setUp()
{
  (*this)->initialize();
}

void 
PyPF::Elemset::setConnectivity(int nelem, int nel, int icone[]) 
{
  if ((*this)->icone == NULL) {
    Elemset::Base* elemset = *this;
    elemset->icone      = new int[nelem*nel];
    if (elemset->elem_conne) delete[] elemset->elem_conne;
    elemset->elem_conne = new int[nel];
  } 
  else {
    Elemset::Base* elemset = *this;
    if (elemset->nelem*elemset->nel != nelem*nel) {
      delete[] elemset->icone; elemset->icone = new int[nelem*nel];
    }
    if (elemset->nel != nel || !elemset->elem_conne) {
      if (elemset->elem_conne) delete[] elemset->elem_conne;
      elemset->elem_conne = new int[nel];
    }
  }

  Elemset::Base* elemset = *this;
  elemset->nelem = nelem;
  elemset->nel   = nel;
  memcpy(elemset->icone, icone, nelem*nel*sizeof(int));
  memset(elemset->elem_conne, 0, nel*sizeof(int));
  elemset->ndof = 0;
}

void 
PyPF::Elemset::getConnectivity(int* nelem, int* nel, int* icone[]) 
{
  if ((*this)->icone == NULL) {
    throw Error("connectivity not set");
  }
  Elemset::Base* elemset = *this;
  *nelem = elemset->nelem;
  *nel   = elemset->nel;
  *icone = elemset->icone;
}

void 
PyPF::Elemset::getSize(int* nelem, int* nel)
{
  if ((*this)->icone == NULL) {
    throw Error("connectivity not set");
  }
  Elemset::Base* elemset = *this;
  *nelem = elemset->nelem;
  *nel   = elemset->nel;
}

int 
PyPF::Elemset::getNDof() 
{
  return (*this)->ndof;
}

void 
PyPF::Elemset::setNDof(int ndof)
{
  (*this)->ndof = ndof;
}


void
PyPF::Elemset::view()
{
  (*this)->print();
}
