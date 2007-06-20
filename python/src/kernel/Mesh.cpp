// $Id$

#include "Mesh.h"
#include "Dofset.h"
#include "Domain.h"

#include <fem.h>
#include <dofmap.h>
#include <readmesh.h>

PF4PY_NAMESPACE_BEGIN
static inline void 
sync(Mesh::Impl* mesh, DTable<double>& nodetable, int ndim)
{
  if (mesh == NULL) return;
  const std::pair<int,int>& shape = nodetable.getShape();
  if (mesh->nodedata == NULL) mesh->nodedata = new ::Nodedata;
  mesh->nodedata->ndim     = ndim;
  mesh->nodedata->nnod     = shape.first;
  mesh->nodedata->nu       = shape.second;
  mesh->nodedata->nodedata = nodetable;
}
typedef pair<int,double*>      FieldEntry;
typedef map<string,FieldEntry> FieldMap;
static inline void 
sync(Mesh::Impl* mesh, const std::string& name, DTable<double>& data)
{
  if (mesh == NULL) return;
  if (mesh->nodedata == NULL) mesh->nodedata = new ::Nodedata;
  FieldMap&   fields = mesh->nodedata->fields;
  FieldEntry& entry = fields[name];
  const std::pair<int,int>& shape = data.getShape();
  entry.first  = shape.second;
  entry.second = data;
}
typedef RefMap<std::string,DTable<double> > MeshFieldMap;
static inline void 
sync(Mesh::Impl* mesh, const MeshFieldMap& fields)
{
  MeshFieldMap::const_iterator i = fields.begin();
  MeshFieldMap::const_iterator e = fields.end();
  while (i != e) { sync(mesh, i->first, *(i->second)); i++; }
}

static inline void
sync(Mesh::Impl* mesh, Elemset::Impl* elemset)
{
  if (mesh == NULL) return;
  if (elemset == NULL) return;
  da_append(mesh->elemsetlist, &elemset);
}
static inline void 
sync(Mesh::Impl* mesh, Options& options)
{
  if (mesh == NULL) return;
  mesh->global_options     = options;
  if (mesh->nodedata != NULL)
    mesh->nodedata->options  = options;
}
PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
class Mesh::Proxy
{
private:
  Proxy();
  Proxy(const Proxy&);
  Proxy& operator=(const Proxy&);
protected:
  std::auto_ptr< Mesh::Impl > _ptr;
public:
  inline operator Mesh::Impl*() const { return this->_ptr.get(); }

public:
  ~Proxy() 
  { 
    Mesh::Impl* mesh = *this;
    if (mesh == NULL) return;
    // + nodedata
    if (mesh->nodedata) {
      mesh->nodedata->nodedata = NULL;
      mesh->nodedata->options  = NULL;
      mesh->nodedata->fields.clear();
    }
    PF4PY_DELETE_SCLR(mesh->nodedata);
    // + elemsetlist
    PF4PY_DELETE_FUNC(da_destroy, mesh->elemsetlist);
    // + options
    mesh->global_options = NULL;
  }

  Proxy(Mesh* M)
    : _ptr(new Mesh::Impl)
  {
    Mesh::Impl* mesh = *this;
    // init mesh members
    mesh->nodedata = NULL;
    mesh->elemsetlist = NULL;
    mesh->global_options = NULL;
    // extract M members
    DTable<double>& nodetable = M->getNodedata();
    int             ndim      = M->ndim;
    int             nnod      = M->nnod;
    int             nu        = nodetable.getShape().second;
    Options&        options   = M->getOptions();
    // create nodedata
    mesh->nodedata = new ::Nodedata;
    // init nodedata members
    mesh->nodedata->ndim     = ndim;
    mesh->nodedata->nnod     = nnod;
    mesh->nodedata->nu       = nu;
    mesh->nodedata->nodedata = nodetable;
    mesh->nodedata->options  = options;
    sync(mesh, M->fields);
    // create elemsetlist
    std::size_t i, n = M->getSize();
    mesh->elemsetlist = da_create_len(sizeof(Elemset::Impl*), n);
    // fill elemsetlist
    if(mesh->elemsetlist == NULL)
      throw Error("Mesh: allocating darray for elemset list");
    for (i=0; i<n; i++) {
      Elemset::Impl* elemset = NULL;
      da_set(mesh->elemsetlist, i, &elemset);
    }
    for (i=0; i<n; i++) {
      Elemset::Impl* elemset = M->getElemset((int)i);
      da_set(mesh->elemsetlist, i, &elemset);
    }
    // set mesh options
    mesh->global_options = mesh->nodedata->options;
  }
}; // class Mesh::Proxy

Mesh::Impl* Mesh::getimpl() const
{ 
  Mesh::Proxy* proxy = this->proxy.get();
  if (proxy != NULL) return *proxy;
  return NULL;
}

PF4PY_NAMESPACE_END


PF4PY_NAMESPACE_BEGIN

Mesh::~Mesh()
{ }

Mesh::Mesh()
  : ndim(1),
    nnod(0),
    nodedata(),
    nodepart(),
    fields(),
    elemsets(),
    options(),
    proxy()
{ 
  throw Error("Mesh::Mesh()"); 
}

Mesh::Mesh(const Mesh& mesh)
  : ndim(mesh.ndim),
    nnod(mesh.nnod),
    nodedata(mesh.nodedata),
    nodepart(mesh.nodepart),
    fields(mesh.fields),
    elemsets(mesh.elemsets),
    options(mesh.options),
    proxy()
{
  throw Error("Mesh::Mesh(const Mesh&)"); 
}

Mesh::Mesh(MPI_Comm comm, int ndim, int nnod)
  : ndim(ndim),
    nnod(nnod),
    nodedata(),
    nodepart((nnod<0)?0:nnod, 1),
    fields(),
    elemsets(),
    options(new Options()),
    proxy()
{
  if (ndim < 1) throw Error("Mesh: ndim < 1");
  if (ndim > 3) throw Error("Mesh: ndim > 3");
  if (nnod < 0) throw Error("Mesh: nnod < 0");
}

DTable<double>& 
Mesh::getNodedata() const
{
  if (!this->nodedata) 
    throw Error("Mesh: nodedata not set");
  return this->nodedata;
}

void
Mesh::setNodedata(const DTable<double>& nodetable)
{
  // check node table
  const std::pair<int,int>& shape = nodetable.getShape();
  if (this->nnod != shape.first)
    throw Error("Mesh: rows != nnod in nodedata");
  if (this->ndim > shape.second)
    throw Error("Mesh: cols < ndim in nodedata");
  // set node table
  this->nodedata = nodetable;
  sync(*this, this->nodedata, this->ndim);
}

DTable<double>& 
Mesh::getField(const std::string& name) const
{
  DTable<double>* entry = this->fields.get(name);
  if (!entry) throw Error("Mesh: field not found");
  return  *entry;
}

void
Mesh::setField(const std::string& name, DTable<double>& data)
{
  const std::pair<int,int>& shape = data.getShape();
  if (shape.first != this->nnod)
    throw Error("Mesh: rows != nnod in data");
  this->fields.set(name, &data);
  sync(*this, name, data);
}

void 
Mesh::getPartitioning(std::vector<int>& part) const
{
  part = this->nodepart;
  for (std::size_t i=0; i<part.size(); part[i++] -= 1);
}

int
Mesh::getSize() const 
{
  return (int) this->elemsets.size();
}

Elemset& 
Mesh::getElemset(int i) const
{
  Elemset* e = this->elemsets.get(static_cast<std::size_t>(i));
  if (e == NULL) 
    throw Error("Mesh: elemset not found, index out of range");
  return *e;
}

void
Mesh::addElemset(const Elemset& elemset)
{
  this->elemsets.add(const_cast<Elemset*>(&elemset));
  sync(*this, elemset);
}

// void
// Mesh::setElemset(int i, Elemset& elemset) 
// {
//   this->elemsets.set(i, &elemset);
//   sync(*this, i, elemset);
// }

// void
// Mesh::delElemset(int i)
// {
//   bool found = this->elemsets.del(i);
//   if (!found) throw Error("elemset not found, index out of range");
// }

Options&
Mesh::getOptions() const
{
  if (!this->options) throw Error("Mesh: options not set");
  return this->options;
}

void
Mesh::setOptions(const Options& options)
{
  this->options = options;
  sync(*this, this->options);
}

PF4PY_NAMESPACE_END

#include "gvars.h"

PF4PY_NAMESPACE_BEGIN

void setup_mesh(MPI_Comm, Mesh::Impl*, std::vector<int>&, std::vector<float>&);

void 
Mesh::setup(Domain* domain)
{
  if (!this->nodedata) 
    throw Error("Mesh: nodedata not set");
  if (this->elemsets.empty())
    throw Error("Mesh: empty elemset list");
  if (!this->options)
    throw Error("Mesh: options not set");

  // create mesh
  this->proxy.reset(new Mesh::Proxy(this));
  Mesh::Impl* mesh = *this;
  int ndof = domain->getNDof();
  // initialize elemsets
  std::size_t i, n = da_length(mesh->elemsetlist);
  for (i=0; i<n; i++) {
    Elemset::Impl** eptr = (Elemset::Impl**) da_ref(mesh->elemsetlist, i);
    assert (eptr != NULL);
    Elemset::Impl* elemset = *eptr;
    assert (elemset != NULL);
    elemset->ndof = ndof;
    elemset->initialize();
  }
  // mesh partitioning
  MPI_Comm            comm  = domain->getComm();
  Options::Impl*      thash = domain->getOptions();
  std::vector<int>&   npart = this->nodepart;
  std::vector<float>& tpwgts= domain->dofset->tpwgts;

  GlobalVars gvars(comm, thash, mesh, NULL);
  setup_mesh(comm, mesh, npart, tpwgts);
}

PF4PY_NAMESPACE_END
