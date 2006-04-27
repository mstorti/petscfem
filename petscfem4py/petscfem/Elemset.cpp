// $Id: Elemset.cpp,v 1.1.2.8 2006/04/27 19:09:17 rodrigop Exp $

#include "Elemset.h"

#include <fem.h>
#include <readmesh.h>

#if 0
void  py_bless_elemset(char*, Elemset*&);
#define bless_elemset py_bless_elemset
#endif

PYPF_NAMESPACE_BEGIN
Elemset::Base*
NewElemset(const std::string& type) {
  Elemset::Base* elemset = NULL;
  bless_elemset(const_cast<char*>(type.c_str()), elemset);
  return elemset;
}
PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN

Elemset::~Elemset()
{ 
  Elemset& handle = *this;
  typedef std::map<std::string,Elemset::Base*> Table;
  Table& table = Elemset::Base::elemset_table;
  Table::iterator m = table.find(handle->name_m);
  if (m != table.end()) table.erase(m);

  Elemset::Base* elemset = *this;
  PYPF_DELETE_VCTR(elemset->type);
  elemset->icone = NULL;
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
  : Handle(NULL, false),
    Object(),
    nelem(0), nel(0), icone()
{ }

Elemset::Elemset(const Elemset& elset)
  : Handle(NULL, false),
    Object(elset),
    nelem(elset.nelem), nel(elset.nel), icone(elset.icone)
{ 
  
  Elemset::Base* input = elset;

  // creation
  std::string type = input->type;
  Elemset::Base* elemset = NewElemset(type.c_str());
  if (elemset == NULL)
    throw Error("unknown elemset type: '" + type + "'");
  elemset->type = new char[type.size()+1];
  strcpy(elemset->type, type.c_str());
  const std::string& name = (Elemset::Base::anon);
  elemset->register_name(name, elemset->type);

  Elemset::Base*& handle = *this;
  handle = elemset;

  // options
  elemset->thash = this->options;

  // connectivity
  elemset->nelem = this->nelem;
  elemset->nel   = this->nel;
  elemset->icone = &this->icone[0];

  elemset->elem_conne = new int[this->nel];
  memset(elemset->elem_conne, 0, this->nel*sizeof(int));

  elemset->ndof = input->ndof;


  int nelem = this->nelem;
  int nel   = this->nel;

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
    elemset->local_store = new void*[input->nelem_here];

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

}

Elemset::Elemset(const std::string& type,
		 int nelem, int nel, const int icone[])
  : Handle(NULL, false), 
    Object(),
    nelem(0), nel(0), icone()
{
  // create
  Elemset::Base* elemset = NewElemset(type.c_str());
  if (elemset == NULL) 
    throw Error("unknown elemset type: '" + type + "'" );
  elemset->type = new char[type.size()+1];
  strcpy(elemset->type, type.c_str());
  const std::string& name = (Elemset::Base::anon);
  elemset->register_name(name, elemset->type);

  Elemset::Base*& handle = *this;
  handle = elemset;

  // options
  elemset->thash = this->options;

  // connectivity
  elemset->nelem = 0;
  elemset->nel   = 0;
  elemset->icone = NULL;
  elemset->ndof  = 0;
  elemset->elem_conne = NULL;

  // partitioning
  elemset->epart       = NULL;
  elemset->epart2      = NULL;
  elemset->e1          = 0;
  elemset->e2          = 0;
  elemset->isfat       = 1;
  elemset->nelem_here  = 0;
  elemset->ghost_elems = da_create(sizeof(int));
  elemset->local_store = NULL;

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

  this->setData(nelem, nel, icone);

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
Elemset::setName(const std::string& name)
{
  PYPF_ASSERT(name.size()>0, "empty name");
  Elemset& elemset = *this;
  typedef std::map<std::string,Elemset::Base*> Table;

  Table& table = Elemset::Base::elemset_table;
  Table::iterator m; 

  m = table.find(elemset->name_m);
  if (m != table.end()) table.erase(m);

  m = table.find(name);
  if (m != table.end()) throw Error("name already registered");
  table[name]     = elemset;
  elemset->name_m = name;
}

void 
Elemset::getDataSize(int* nelem, int* nel) const
{
  Elemset::Base* elemset = *this;
  if (nelem) *nelem = elemset->nelem;
  if (nel)   *nel   = elemset->nel;
}

void 
Elemset::getData(int* nelem, int* nel, const int* icone[]) const
{
  Elemset::Base* elemset = *this;
  if (nelem) *nelem = this->nelem;
  if (nel)   *nel   = this->nel;
  if (icone) *icone = &this->icone[0];
}

#undef  PYPF_NEW_ARRAY
#define PYPF_NEW_ARRAY(array, old_size, new_size)       \
do {                                                    \
    if ((array) == NULL && (new_size)>0) {              \
      (array) = new int[(new_size)];                    \
    } else if ((old_size) != (new_size)) {              \
      delete[] (array); (array) = new int[(new_size)];  \
    }                                                   \
} while(0)

void 
Elemset::setData(int nelem, int nel, int const icone[])
{
  // test data
  if (nelem*nel != 0) {
    PYPF_ASSERT(nelem>=1,    "invalid number of nodes (nelem<1)");
    PYPF_ASSERT(nel>=1,      "invalid number of element nodes (nel<1)");
    PYPF_ASSERT(icone!=NULL, "null pointer to data array");
  } else nelem = nel = 0;
  if (this->nelem*this->nel != 0) {
    PYPF_ASSERT(this->nelem==nelem, "cannot change data size (nelem)");
    PYPF_ASSERT(this->nel==nel,     "cannot change data size (nel)");
  }

  this->nelem = nelem;
  this->nel   = nel;
  this->icone.resize(nelem*nel);
  memcpy(&this->icone[0], icone, nelem*nel*sizeof(int));
  for (int i=0; i<nelem*nel; this->icone[i++]+=1);


  Elemset::Base* elemset = *this;
  elemset->nelem = this->nelem;
  elemset->nel   = this->nel;
  elemset->icone = &this->icone[0];
  PYPF_NEW_ARRAY(elemset->elem_conne, elemset->nel,                nel);
  PYPF_NEW_ARRAY(elemset->epart,      elemset->nelem,              nelem);
  PYPF_NEW_ARRAY(elemset->epart2,     elemset->nelem,              nelem);
  memset(elemset->elem_conne, 0, sizeof(int)*nel);
  memset(elemset->epart,      0, sizeof(int)*nelem);
  memset(elemset->epart2,     0, sizeof(int)*nelem);

  PYPF_DELETE_FUNC(da_destroy, elemset->ghost_elems);
  elemset->ghost_elems = da_create(sizeof(int));
}

#undef PYPF_NEW_ARRAY


void
Elemset::getPart(int* n, int* part[]) const
{
  Elemset::Base* elemset = *this;
  if (n)    *n    = elemset->nelem;
  if (part) *part = elemset->epart;
}


void
Elemset::getElem(int i, int* n, const int* elem[]) const
{
  int nelem, nel; const int* icone;
  this->getData(&nelem, &nel, &icone);
  PYPF_ASSERT(i>=0 && i<nelem, "index out of range");
  if (n)    *n    = nel;
  if (elem) *elem = &icone[i*nel];
}

void 
Elemset::setElem(int i, int n, const int elem[])
{
  int nelem, nel; const int* icone;
  this->getData(&nelem, &nel, &icone);
  PYPF_ASSERT(i>=0 && i<nelem, "index out of range");
  PYPF_ASSERT(n != nel, "invalid data size");
  int* data = const_cast<int*>(&icone[i*nel]);
  for (int k=0; k<nel; k++) data[k] = elem[k]+1;
}

int Elemset::getSize() const
{
  int nelem;
  this->getDataSize(&nelem, NULL);
  return nelem;
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
  PYPF_ASSERT(ndof>=1, "value out of range, (ndof<1)");
  Elemset::Base* elemset = *this;
  elemset->ndof = ndof;
}

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
