// $Id$

#ifndef PF4PY_MESH_H
#define PF4PY_MESH_H

#include <memory>
#include "petscfem4py.h"

#include "Comm.h"
#include "Object.h"
#include "Options.h"
#include "DTable.h"
#include "Elemset.h"


PF4PY_NAMESPACE_BEGIN

class Mesh 
  : public Object
{
  friend class Dofset;
  friend class Domain;

public:
  typedef ::Mesh Impl;

protected:
  int               ndim;
  int               nnod;
  std::vector<int>  nodepart;

  RefVal<DTable<double> >             nodedata;
  RefMap<std::string,DTable<double> > fields;
  RefVec<Elemset>                     elemsets;
  RefVal<Options>                     options;

#if !defined(SWIG)
protected:
  class Proxy; friend class Proxy;
  std::auto_ptr< Proxy > proxy;
  Mesh::Impl* getimpl() const;
public:
  inline operator Mesh::Impl*() const { return this->getimpl(); }
#endif

protected:
  void setup(Domain*);

public:
  ~Mesh();
private:
  Mesh();
  Mesh(const Mesh& mesh);
  Mesh& operator=(const Mesh& mesh);
protected:
  Mesh(MPI_Comm comm, int ndim, int nnod);

public:
  DTable<double>& getNodedata() const;
  void            setNodedata(const DTable<double>& nodedata);

  DTable<double>& getField(const std::string& name) const;
  void            setField(const std::string& name,
			   DTable<double>& data);

public:
  void getPartitioning(std::vector<int>& part) const;

public:
  int      getSize() const;
  Elemset& getElemset(int index) const;
  void     addElemset(const Elemset& elemset);

public:
  Options& getOptions() const;
  void     setOptions(const Options& options);
  
  
};

PF4PY_NAMESPACE_END

#endif // PF4PY_MESH_H

// Local Variables:
// mode: C++
// End:
