// $Id$

#ifndef PF4PY_ELEMSET_H
#define PF4PY_ELEMSET_H

#include <memory>
#include <string>
#include <vector>
#include <map>
#include "petscfem4py.h"

#include "Object.h"
#include "DTable.h"
#include "PTable.h"
#include "Options.h"

PF4PY_NAMESPACE_BEGIN

class Elemset
  : public Object
{
  friend class Mesh;
  friend class Domain;

public:
  typedef ::Elemset Impl;
  
protected:
  std::string             type;
  RefVal<DTable<int> >    conntable;
  RefVal<PTable<int> >    proptable_i;
  RefVal<PTable<double> > proptable_s;
  RefVal<Options>         options;

#if !defined(SWIG)
protected:
  class Proxy; friend class Proxy;
  std::auto_ptr< Proxy > proxy;
  Elemset::Impl* getimpl() const;
public:
  inline operator Elemset::Impl*() const { return this->getimpl(); }
#endif

protected:
  void setup(Domain*);

public:
  ~Elemset();
private:
  Elemset(const Elemset&);
  Elemset& operator=(const Elemset&);
public:
  Elemset();
  Elemset(const std::string& type);
  Elemset(const std::string& type, const DTable<int>& icone);
  Elemset(const std::string& type, const DTable<int>& icone, const Options& options);

public:
  const std::string& getType() const;
  void setType(const std::string& type);

  Options& getOptions() const;
  void setOptions(const Options& options);

  DTable<int>& getData() const;
  void setData(const DTable<int>& conntable);
  
  void setPTable(const PTable<int>&    proptable);
  PTable<int>& getPTableI() const
  { 
    if (!this->proptable_i) throw Error("Elemset: integer property table not set");
    return this->proptable_i;
  }

  void setPTable(const PTable<double>& proptable);
  PTable<double>& getPTableS() const
  { 
    if (!this->proptable_s) throw Error("Elemset: scalar property table not set");
    return this->proptable_s;
  }
  
};


PF4PY_NAMESPACE_END

#endif // PF4PY_ELEMSET_H

// Local Variables:
// mode: C++
// End:
