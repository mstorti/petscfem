// $Id: Elemset.cpp,v 1.1.2.5 2006/03/20 16:06:00 rodrigop Exp $

#include "Elemset.h"

#include <fem.h>
#include <readmesh.h>

PYPF_NAMESPACE_BEGIN

OptionTable*
Elemset::get_opt_table() const
{ 
  Elemset::Base* elemset = *this;
  if (elemset->thash == NULL) elemset->thash = new OptionTable;
  return elemset->thash;
}

Elemset::~Elemset()
{ 
  Elemset::Base* elemset = *this;
  PYPF_DELETE_VCTR(elemset->type);
  PYPF_DELETE_VCTR(elemset->icone);
  PYPF_DELETE_VCTR(elemset->epart);
  PYPF_DELETE_VCTR(elemset->epart2);
  PYPF_DELETE_VCTR(elemset->elemprops);
  PYPF_DELETE_VCTR(elemset->elemiprops);
  PYPF_DELETE_VCTR(elemset->elemprops_add);
  PYPF_DELETE_VCTR(elemset->elemiprops_add);
  PYPF_DELETE_SCLR(elemset->thash);
  PYPF_DELETE(g_hash_table_destroy,elemset->elem_iprop_names);
  PYPF_DELETE(g_hash_table_destroy,elemset->elem_prop_names);
  PYPF_DELETE_VCTR(elemset->elem_conne);
  PYPF_DELETE_SCLR(elemset);
}

Elemset::Elemset()
  : Ptr(NULL, false), Object()
{ 
#if 0
  Elemset::Base* elemset = *this;
  elemset->type             = NULL;
  elemset->icone            = NULL;
  elemset->epart            = NULL;
  elemset->epart2           = NULL;
  elemset->elemprops        = NULL;
  elemset->elemiprops       = NULL;
  elemset->elemprops_add    = NULL;
  elemset->elemiprops_add   = NULL;
  elemset->elem_iprop_names = NULL;
  elemset->elem_prop_names  = NULL;
  elemset->elem_conne       = NULL;
  elemset->thash            = NULL;
#endif
}

Elemset::Elemset(const Elemset& _elemset)
  : Ptr(NULL, false), Object(_elemset)
{ 
  
  Elemset::Base* input = _elemset;

  // create
  Elemset::Base* elemset = NULL;
  char* type = new char[strlen(input->type)+1];
  strcpy(type, input->type);
  bless_elemset(type, elemset);
  if (elemset == NULL) 
    { delete[] type; throw Error("unknown elemset type"); }
  elemset->type = type;
  const std::string& name = (Elemset::Base::anon);
  elemset->register_name(name, elemset->type);
  
  this->ptr = elemset;

  // connectivity
  int nelem, nel; int* icone;
  _elemset.getConnectivity(&nelem, &nel, &icone);
  this->setConnectivity(nelem, nel, icone);

  // partitioning
  elemset->ndof  = input->ndof;
  elemset->isfat = input->isfat;
  memcpy(elemset->epart, input->epart, elemset->nelem * sizeof(int));

  // options
  elemset->thash = new OptionTable;
  this->setOptions(_elemset.getOptions());
}


Elemset::Elemset(Elemset::Base* _elemset)
  : Ptr(_elemset), Object()
{ }

Elemset::Elemset(const std::string& type,
		 const std::string& name)
  : Ptr(NULL, false), Object()
{
  // create
  Elemset::Base* elemset = NULL;
  char* _type = new char[type.size()+1];
  strcpy(_type, type.c_str());
  bless_elemset(_type, elemset);
  if (elemset == NULL) {
    delete[] _type; throw Error("unknown elemset type");
  }
  elemset->type = _type;
  const std::string& _name =  (name.size()) ? name : (Elemset::Base::anon);
  elemset->register_name(_name, elemset->type);

  this->ptr = elemset;

  elemset->ndof  = 0;
  elemset->isfat = 0;
  
  //  connectivity
  elemset->nelem      = 0;
  elemset->nel        = 0;
  elemset->icone      = NULL;
  elemset->elem_conne = NULL;

  // partitioning
  elemset->epart = NULL;

  // options
  elemset->thash = new OptionTable;

  // int props
  elemset->neliprops        = 0; 
  elemset->elemiprops       = NULL; 
  elemset->elem_iprop_names = g_hash_table_new(&g_str_hash,&g_str_equal); 
  elemset->neliprops_add    = 0; 
  elemset->elemiprops_add   = NULL; 

  // double props
  elemset->nelprops         = 0; 
  elemset->elemprops        = NULL;
  elemset->elem_prop_names  = g_hash_table_new(&g_str_hash,&g_str_equal); 
  elemset->nelprops_add     = 0; 
  elemset->elemprops_add    = NULL; 
}

std::string 
Elemset::getType() const
{
  Elemset::Base* elemset = *this;
  return elemset->type;
}

std::string 
Elemset::getName() const
{
  Elemset::Base* elemset = *this;
  return elemset->name();
}

void 
Elemset::getConnectivity(int* nelem, int* nel, int* icone[]) const
{
  Elemset::Base* elemset = *this;
  *nelem = elemset->nelem;
  *nel   = elemset->nel;
  *icone = elemset->icone;
}

void 
Elemset::setConnectivity(int nelem, int nel, int icone[]) 
{

  Elemset::Base* elemset = *this;

  if (elemset->icone == NULL) {
    elemset->icone = new int[nelem*nel];
  }
  else if (elemset->nelem*elemset->nel != nelem*nel) {
    delete[] elemset->icone; 
    elemset->icone = new int[nelem*nel];
  }
  if (elemset->elem_conne == NULL) {
    elemset->elem_conne = new int[nel];
  } 
  else if (elemset->nel != nel) {
    delete[] elemset->elem_conne;
    elemset->elem_conne = new int[nel];
  }
  if (elemset->epart == NULL) {
    elemset->epart = new int[nelem];
  } else if (elemset->nelem != nelem) {
    delete[] elemset->epart;
    elemset->epart = new int[nelem];
  }
  
  elemset->nelem = nelem;
  elemset->nel   = nel;
  memcpy(elemset->icone,      icone, nelem*nel * sizeof(int));
  memset(elemset->elem_conne,     0, nel       * sizeof(int));
  memset(elemset->epart,          0, nelem     * sizeof(int));
}

void 
Elemset::getSize(int* nelem, int* nel) const
{
  Elemset::Base* elemset = *this;
  *nelem = elemset->nelem;
  *nel   = elemset->nel;
}

int 
Elemset::getNDof() const
{
  Elemset::Base* elemset = *this;
  return elemset->ndof;
}

void 
Elemset::setNDof(int ndof)
{
  Elemset::Base* elemset = *this;
  elemset->ndof = ndof;
}

void 
Elemset::setUp()
{
  Elemset::Base* elemset = *this;
  elemset->initialize();
}

void
Elemset::clear()
{
  Elemset::Base* elemset = *this;
  //PYPF_DELETE_VCTR(elemset->type);
  PYPF_DELETE_VCTR(elemset->icone);
  PYPF_DELETE_VCTR(elemset->epart);
  PYPF_DELETE_VCTR(elemset->epart2);
  PYPF_DELETE_VCTR(elemset->elemprops);
  PYPF_DELETE_VCTR(elemset->elemiprops);
  PYPF_DELETE_VCTR(elemset->elemprops_add);
  PYPF_DELETE_VCTR(elemset->elemiprops_add);
  PYPF_DELETE_SCLR(elemset->thash);
  PYPF_DELETE(g_hash_table_destroy,elemset->elem_iprop_names);
  PYPF_DELETE(g_hash_table_destroy,elemset->elem_prop_names);
  PYPF_DELETE_VCTR(elemset->elem_conne);
  //PYPF_DELETE_SCLR(elemset);
}

void
Elemset::view() const
{
  Elemset::Base* elemset = *this;
  elemset->print();
}

PYPF_NAMESPACE_END
