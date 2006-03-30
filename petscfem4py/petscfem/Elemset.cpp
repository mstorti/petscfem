// $Id: Elemset.cpp,v 1.1.2.7 2006/03/30 15:18:14 rodrigop Exp $

#include "Elemset.h"

#include <fem.h>
#include <readmesh.h>

PYPF_NAMESPACE_BEGIN

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
  elemset->thash = NULL;
  PYPF_DELETE_FUNC(g_hash_table_destroy,elemset->elem_prop_names);
  PYPF_DELETE_FUNC(g_hash_table_destroy,elemset->elem_iprop_names);
  PYPF_DELETE_FUNC(da_destroy, elemset->ghost_elems);
  PYPF_DELETE_VCTR(elemset->local_store);
  PYPF_DELETE_VCTR(elemset->elem_conne);
  delete elemset;
}

Elemset::Elemset()
  : Handle(NULL, false), Object(),
    nelem(0), nel(0), icone()
{ }

Elemset::Elemset(const Elemset& elset)
  : Handle(NULL, false), Object(elset),
    nelem(elset.nelem), nel(elset.nel), icone(elset.icone)
{ 
  
  Elemset::Base* input = elset;

  // creation
  Elemset::Base* elemset = NULL;
  bless_elemset(input->type, elemset);
  if (elemset == NULL) throw Error("unknown elemset type");
  elemset->type = new char[strlen(input->type)+1];
  strcpy(elemset->type, input->type);
  const std::string& name = (Elemset::Base::anon);
  elemset->register_name(name, elemset->type);
  
  Elemset::Base*& handle = *this;
  handle = elemset;

  // connectivity
  int nelem = elemset->nelem = input->nelem;
  int nel   = elemset->nel   = input->nel;
  elemset->icone  = new int[nelem*nel];
  memcpy(elemset->icone,  input->icone,  sizeof(int)*nelem*nel);
  elemset->elem_conne = new int[nel];
  memset(elemset->elem_conne, 0, nel*sizeof(int));

  elemset->ndof = input->ndof;

  // partitioning
  elemset->epart       = new int[nelem];
  memcpy(elemset->epart, input->epart,  sizeof(int)*nelem);
  elemset->epart2      = new int[nelem];
  memcpy(elemset->epart2, input->epart2, sizeof(int)*nelem);
  elemset->epart_p     = input->epart_p;
  elemset->e1          = input->e1;
  elemset->e2          = input->e2;
  elemset->isfat       = input->isfat;
  elemset->nelem_here  = input->nelem_here;
  elemset->ghost_elems = da_create(sizeof(int));
  da_concat_da(elemset->ghost_elems, input->ghost_elems);
  elemset->local_store = NULL;
  if (input->local_store)
    elemset->local_store = new (void*)[input->nelem_here];

  // options
  elemset->thash = this->options;

  // double properties
  int nelprops = elemset->nelprops = input->nelprops;
  elemset->elemprops = new double[nelem*nelprops];
  memcpy(elemset->elemprops, input->elemprops, sizeof(double)*nelem*nelprops);
  elemset->elem_prop_names  = g_hash_table_new(&g_str_hash,&g_str_equal);
  
  int nelprops_add = elemset->nelprops_add = input->nelprops_add;
  elemset->elemprops_add = new double[elemset->nelem*nelprops_add];
  memcpy(elemset->elemprops_add, input->elemprops_add, sizeof(double)*nelem*nelprops_add);

  // int properties
  int neliprops = elemset->neliprops = input->neliprops;
  elemset->elemiprops = new int[nelem*neliprops];
  memcpy(elemset->elemiprops, input->elemiprops, sizeof(int)*nelem*neliprops);
  elemset->elem_iprop_names = g_hash_table_new(&g_str_hash,&g_str_equal);

  int neliprops_add = elemset->neliprops_add = input->neliprops_add;
  elemset->elemiprops_add = new int[elemset->nelem*neliprops_add];
  memcpy(elemset->elemiprops_add, input->elemiprops_add, sizeof(int)*nelem*neliprops_add);

  
  elemset->initialize();
  
}


Elemset::Elemset(Elemset::Base* base)
  : Handle(base), Object(),
    nelem(base->nelem), nel(base->nel)
{ 
  Elemset::Base* elemset = *this;
  // options
  if (elemset->thash == NULL) 
    elemset->thash = this->options;
  else 
    this->options = elemset->thash;
}

Elemset::Elemset(const std::string& type,
		 const std::string& name)
  : Handle(NULL, false), Object(),
    nelem(0), nel(0), icone()
{
  // create
  Elemset::Base* elemset = NULL;
  bless_elemset(const_cast<char*>(type.c_str()), elemset);
  if (elemset == NULL) throw Error("unknown elemset type");
  elemset->type = new char[type.size()+1];
  strcpy(elemset->type, type.c_str());
  const std::string& _name =  (name.size()) ? name : (Elemset::Base::anon);
  elemset->register_name(_name, elemset->type);

  Elemset::Base*& handle = *this;
  handle = elemset;

  //  connectivity
  elemset->nelem      = 0;
  elemset->nel        = 0;
  elemset->icone      = NULL;
  elemset->elem_conne = NULL;

  elemset->ndof  = 0;
  
  // partitioning
  elemset->epart       = NULL;
  elemset->epart2      = NULL;
  elemset->e1          = 0;
  elemset->e2          = 0;
  elemset->isfat       = 1;
  elemset->nelem_here  = 0;
  elemset->ghost_elems = da_create(sizeof(int));
  elemset->local_store = NULL;

  // options
  elemset->thash = this->options;

  // double props
  elemset->nelprops         = 0;
  elemset->elemprops        = NULL;
  elemset->elem_prop_names  = g_hash_table_new(&g_str_hash,&g_str_equal); 
  elemset->nelprops_add     = 0;
  elemset->elemprops_add    = NULL;

  // int props
  elemset->neliprops        = 0; 
  elemset->elemiprops       = NULL; 
  elemset->elem_iprop_names = g_hash_table_new(&g_str_hash,&g_str_equal); 
  elemset->neliprops_add    = 0; 
  elemset->elemiprops_add   = NULL; 

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
Elemset::getSize(int* nelem, int* nel) const
{
  Elemset::Base* elemset = *this;
  if (nelem) *nelem = elemset->nelem;
  if (nel)   *nel   = elemset->nel;
}

void 
Elemset::getConnectivity(int* nelem, int* nel, int* icone[]) const
{
  Elemset::Base* elemset = *this;
  if (nelem) *nelem = elemset->nelem;
  if (nel)   *nel   = elemset->nel;
  if (icone) *icone = elemset->icone;
}

#undef  PYPF_NEW_ARRAY
#define PYPF_NEW_ARRAY(array, old_size, new_size)       \
do {                                                    \
    if ((array) == NULL) {                              \
      (array) = new int[(new_size)];                    \
    } else if ((old_size) != (new_size)) {              \
      delete[] (array); (array) = new int[(new_size)];  \
    }                                                   \
} while(0)

void 
Elemset::setConnectivity(int nelem, int nel, int icone[])
{
  Elemset::Base* elemset = *this;

  PYPF_NEW_ARRAY(elemset->icone,      elemset->nelem*elemset->nel, nelem*nel);
  PYPF_NEW_ARRAY(elemset->elem_conne, elemset->nel,                nel);
  PYPF_NEW_ARRAY(elemset->epart,      elemset->nelem,              nelem);
  PYPF_NEW_ARRAY(elemset->epart2,     elemset->nelem,              nelem);

  this->nelem = elemset->nelem = nelem;
  this->nel   = elemset->nel   = nel;
  memcpy(elemset->icone, icone, sizeof(int)*nelem*nel);
  for (int i=0; i<nelem*nel; elemset->icone[i++]+=1);
  memset(elemset->elem_conne, 0, sizeof(int)*nel);
  memset(elemset->epart,      0, sizeof(int)*nelem);
  memset(elemset->epart2,     0, sizeof(int)*nelem);

  PYPF_DELETE_FUNC(da_destroy, elemset->ghost_elems);
  elemset->ghost_elems = da_create(sizeof(int));
}

#undef PYPF_NEW_ARRAY

std::vector<int> 
Elemset::getElem(int n) const 
{
  int nelem, nel; int* icone;
  this->getConnectivity(&nelem, &nel, &icone);
  if (n<0 || n>=nelem) throw Error("index out of range");
  const int* data = &icone[n*nel];
  std::vector<int> elem(nel);
  for (int i=0; i<nel; i++) elem[i] = data[i]-1;
  return elem;
}

void 
Elemset::setElem(int n, const std::vector<int>& elem)
{
  int nelem, nel; int* icone;
  this->getConnectivity(&nelem, &nel, &icone);
  if (n<0 || n>=nelem) throw Error("index out of range");
  if (elem.size() != nel) throw Error("invalid number of dimensions");
  int* data = &icone[n*nel];
  for (int i=0; i<nel; i++) data[i] = elem[i]+1;
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
{ }

void
Elemset::clear()
{
  Elemset::Base* elemset = *this;

  // options
  this->options.clear();
  elemset->thash = this->options;

  // connectivity
  elemset->nelem = 0;
  elemset->nel   = 0;
  PYPF_DELETE_VCTR(elemset->icone);
  PYPF_DELETE_VCTR(elemset->elem_conne);

  elemset->ndof  = 0;
  
  // partition
  PYPF_DELETE_VCTR(elemset->epart);
  PYPF_DELETE_VCTR(elemset->epart2);
  elemset->e1          = 0;
  elemset->e2          = 0;
  elemset->isfat       = 1;
  elemset->nelem_here  = 0;
  PYPF_DELETE_VCTR(elemset->local_store);
  PYPF_DELETE_FUNC(da_destroy, elemset->ghost_elems);
  elemset->ghost_elems = da_create(sizeof(int));
  
  // properties
  elemset->nelprops      = 0;
  elemset->neliprops     = 0;
  elemset->nelprops_add  = 0;
  elemset->neliprops_add = 0;
  PYPF_DELETE_VCTR(elemset->elemprops);
  PYPF_DELETE_VCTR(elemset->elemiprops);
  PYPF_DELETE_VCTR(elemset->elemprops_add);
  PYPF_DELETE_VCTR(elemset->elemiprops_add);
  PYPF_DELETE_FUNC(g_hash_table_destroy, elemset->elem_prop_names);
  PYPF_DELETE_FUNC(g_hash_table_destroy, elemset->elem_iprop_names);
  elemset->elem_prop_names = g_hash_table_new(&g_str_hash,&g_str_equal); 
  elemset->elem_prop_names = g_hash_table_new(&g_str_hash,&g_str_equal); 
}

void
Elemset::view() const
{
  Elemset::Base* elemset = *this;
  elemset->print();
}

PYPF_NAMESPACE_END
