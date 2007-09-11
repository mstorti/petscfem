// $Id$

#include "Elemset.h"
#include "Domain.h"

#include <fem.h>
#include <readmesh.h>

PF4PY_NAMESPACE_BEGIN
static ::Elemset*
destroy(::Elemset* elemset)
{
  if (elemset == NULL) return NULL;

  typedef std::map<std::string,::Elemset*> ElemTable;
  ElemTable& table = ::Elemset::elemset_table;
  ElemTable::iterator m = table.find(elemset->name_m);
  if (m != table.end()) table.erase(m);

  PF4PY_DELETE_VCTR(elemset->type);
  elemset->thash = NULL;
  elemset->icone = NULL;
  PF4PY_DELETE_VCTR(elemset->elem_conne);

  PF4PY_DELETE_VCTR(elemset->epart);
  PF4PY_DELETE_VCTR(elemset->epart2);
  PF4PY_DELETE_FUNC(da_destroy, elemset->ghost_elems);
  PF4PY_DELETE_VCTR(elemset->local_store);

  elemset->elemprops  = NULL;
  elemset->elemiprops = NULL;
  PF4PY_DELETE_FUNC(g_hash_table_destroy, elemset->elem_prop_names);
  PF4PY_DELETE_FUNC(g_hash_table_destroy, elemset->elem_iprop_names);

  PF4PY_DELETE_VCTR(elemset->elemprops_add);
  PF4PY_DELETE_VCTR(elemset->elemiprops_add);

  return elemset;
}
static void
sync(::Elemset* elemset, const std::string& type, const std::string& name)
{
  if (elemset == NULL) return;
  elemset->type = new char[type.size()+1];
  strcpy(elemset->type, type.c_str());
  elemset->register_name(name, elemset->type);
}
static void
sync(::Elemset* elemset, DTable<int>& conntable)
{ 
  if (elemset == NULL) return;
  const std::vector<int>&   array = conntable.getArray();
  const std::pair<int,int>& shape = conntable.getShape();
  int  nelem = shape.first;
  int  nel   = shape.second;
  int* icone = const_cast<int*>(&array[0]);
  elemset->icone = icone;
  elemset->nelem = nelem;
  elemset->nel   = nel;
  elemset->ndof  = 0;

  PF4PY_DELETE_VCTR(elemset->elem_conne);
  PF4PY_DELETE_VCTR(elemset->epart);
  PF4PY_DELETE_VCTR(elemset->epart2);
  elemset->elem_conne = new int[nel];
  elemset->epart      = new int[nelem];
  elemset->epart2     = new int[nelem];
  memset(elemset->elem_conne, 0, sizeof(int)*nel);
  memset(elemset->epart,      0, sizeof(int)*nelem);
  memset(elemset->epart2,     0, sizeof(int)*nelem);

  PF4PY_DELETE_FUNC(da_destroy, elemset->ghost_elems);
  elemset->ghost_elems = da_create(sizeof(int));
  if (elemset->ghost_elems == NULL)
    throw Error("Elemset: allocating darray for ghost elements");
}
static void
sync(::Elemset* elemset, PTable<int>& ptable)
{ 
  if (elemset == NULL) return;
  elemset->neliprops  = 0;
  elemset->elemiprops = NULL;
  PF4PY_DELETE_FUNC(g_hash_table_destroy, elemset->elem_iprop_names);
  elemset->elem_iprop_names = g_hash_table_new(&g_str_hash, &g_str_equal);
  if (elemset->elem_iprop_names == NULL)
    throw Error("Elemset: allocating g_hash_table for integer properties");
}
static void
sync(::Elemset* elemset, PTable<double>& ptable)
{ 
  if (elemset == NULL) return;
  elemset->nelprops  = 0;
  elemset->elemprops = NULL;
  PF4PY_DELETE_FUNC(g_hash_table_destroy, elemset->elem_prop_names);
  elemset->elem_prop_names  = g_hash_table_new(&g_str_hash, &g_str_equal); 
  if (elemset->elem_prop_names == NULL)
    throw Error("Elemset: allocating g_hash_table for scalar properties");
}
static void
sync(::Elemset* elemset, Options& options)
{ 
  if (elemset == NULL) return;
  elemset->thash = options;
}
PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
class Elemset::Proxy
{
private:
  Proxy();
  Proxy(const Proxy&);
  Proxy& operator=(const Proxy&);
protected:
  std::auto_ptr< Elemset::Impl > _ptr;
public:
  inline operator Elemset::Impl*() const { return this->_ptr.get(); }

public:
  ~Proxy() 
  {
    Elemset::Impl* elemset = *this;
    if (elemset == NULL) return;
    destroy(elemset);
  }

  Proxy(Elemset* M, Elemset::Impl* impl)
    : _ptr(impl)
  {
    Elemset::Impl* elemset = *this;
    if (elemset == NULL) return;
  }
}; // class Elemset::Proxy

Elemset::Impl* Elemset::getimpl() const
{ 
  Elemset::Proxy* proxy = this->proxy.get();
  if (proxy != NULL) return *proxy;
  return NULL;
}

PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN

Elemset::~Elemset() 
{  
  destroy(*this); 
}

Elemset::Elemset()
  : type(),
    conntable(),
    proptable_i(),
    proptable_s(),
    options(new Options()),
    proxy()
{
}

Elemset::Elemset(const std::string& type)
  : type(type),
    conntable(),
    proptable_i(),
    proptable_s(),
    options(new Options()),
    proxy()
{
}

Elemset::Elemset(const std::string& type, 
		 const DTable<int>& conntable)
  : type(type),
    conntable(conntable),
    proptable_i(),
    proptable_s(),
    options(new Options()),
    proxy()
{ 
}

Elemset::Elemset(const std::string& type, 
		 const DTable<int>& conntable, 
		 const Options& options)
  : type(type),
    conntable(conntable),
    proptable_i(),
    proptable_s(),
    options(options),
    proxy()
{ 
}

const std::string&
Elemset::getType() const
{
  if (!this->type.size()) throw Error("Elemset: type not set");
  return this->type;
}

void
Elemset::setType(const std::string& type)
{
  if (this->type.size()) throw Error("Elemset: cannot change type");
  this->type = type;
}


DTable<int>&
Elemset::getData() const
{
  if (!this->conntable) throw Error("Elemset: connectivity not set");
  return this->conntable;
}

void
Elemset::setData(const DTable<int>& conntable)
{
  if (!this->conntable)
    this->conntable = conntable;
  else
    throw Error("Elemset: cannot change connectivity");
  sync(*this, this->conntable);
}

void 
Elemset::setPTable(const PTable<int>& pt)
{
  // check
  if (!this->conntable) throw Error("Elemset: connectivity not set");
  if (pt.getShape().first != this->conntable->getShape().first)
    throw Error("Elemset: rows != nelem in integer property table");
  // set property table
  this->proptable_i = pt;
  sync(*this, this->proptable_i);
}

void 
Elemset::setPTable(const PTable<double>& pt)
{
  // check
  if (!this->conntable) throw Error("Elemset: connectivity not set");
  if (pt.getShape().first != this->conntable->getShape().first)
    throw Error("Elemset: rows != nelem in scalar property table");
  // set property table
  this->proptable_s = pt;
  sync(*this, this->proptable_s);
}

Options&
Elemset::getOptions() const
{
  if (!this->options) throw Error("Elemset: options not set");
  return this->options;
}

void
Elemset::setOptions(const Options& options)
{
  this->options = options;
  sync(*this, this->options);
}

PF4PY_NAMESPACE_END


void bless_elemset(char *t,Elemset *& e) { e=NULL; }

void bless_elemset0(char *,Elemset *&);
void bless_elemset_ns(char *,Elemset *&);
void bless_elemset_advdif(char *,Elemset *&);


PF4PY_NAMESPACE_BEGIN

static Elemset::Impl* init(Elemset::Impl* elemset)
{
  if (elemset == NULL) return NULL;
  // type
  elemset->type = NULL;
  // options
  elemset->thash = NULL;
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
  elemset->ghost_elems = NULL;
  elemset->local_store = NULL;
  // scalar props
  elemset->nelprops         = 0;
  elemset->elemprops        = NULL;
  elemset->elem_prop_names  = NULL;
  elemset->nelprops_add     = 0;
  elemset->elemprops_add    = NULL;
  // int props
  elemset->neliprops        = 0; 
  elemset->elemiprops       = NULL; 
  elemset->elem_iprop_names = NULL;
  elemset->neliprops_add    = 0; 
  elemset->elemiprops_add   = NULL; 
  // prop tables
  elemset->elem_prop_names  = NULL; 
  elemset->elem_iprop_names = NULL; 
  //
  return elemset;
}

Elemset::Impl* bless_ns(const std::string& type)
{
  Elemset::Impl* elemset = NULL;
  char *etype = const_cast<char*>(type.c_str());
  bless_elemset_ns(etype, elemset);
  return init(elemset);
}

Elemset::Impl* bless_ad(const std::string& type)
{
  Elemset::Impl* elemset = NULL;
  char *etype = const_cast<char*>(type.c_str());
  bless_elemset_advdif(etype, elemset);
  return init(elemset);
}

void 
Elemset::setup(Domain* domain)
{
  if (!this->type.size()) throw Error("Elemset: type not set");
  if (!this->conntable)   throw Error("Elemset: connectivity not set");
  if (!this->options)     throw Error("Elemset: options not set");
  // create elemset
  Elemset::Impl* elemset = 0;
  const std::string app = domain->getType();
  if      (app == "NS") elemset = bless_ns(this->type);
  else if (app == "AD") elemset = bless_ad(this->type);
  else throw Error("Elemset: invalid domain application type");
  if (elemset == NULL) throw Error("Elemset: unknown type '" + this->type +
				   "' for application '" + app + "'");
  this->proxy.reset(new Elemset::Proxy(this, elemset));
  // initialize elemset
  sync(elemset, this->type, (::Elemset::anon));
  sync(elemset, this->conntable);
  sync(elemset, this->proptable_i);
  sync(elemset, this->proptable_s);
  sync(elemset, this->options);
  elemset->ndof = domain->getNDof();
  elemset->initialize();
}

PF4PY_NAMESPACE_END
