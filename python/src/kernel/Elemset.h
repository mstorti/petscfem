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
  std::auto_ptr< Elemset::Impl > impl;
public:
  inline operator Elemset::Impl*() const { return this->impl.get(); }
#endif

public:
  ~Elemset();
private:
  Elemset();
  Elemset(const Elemset&);
  Elemset& operator=(const Elemset&);
public:
  Elemset(const std::string& type);
  Elemset(const std::string& type, const DTable<int>& icone);
  Elemset(const std::string& type, const DTable<int>& icone, const Options& options);

public:
  const std::string& getType() const;

  Options& getOptions() const;
  void setOptions(const Options& options);

  DTable<int>& getData() const;
  void setData(const DTable<int>& conntable);

  
  void setPTable(const PTable<int>&    proptable);
  void setPTable(const PTable<double>& proptable);
  
  PTable<int>& Elemset::getPTableI() const
  { 
    if (!this->proptable_i) 
      throw Error("Elemset: integer property table not set");
    return this->proptable_i;
  };
  PTable<double>& Elemset::getPTableS() const
  { 
    if (!this->proptable_s) 
      throw Error("Elemset: scalar property table not set");
    return this->proptable_s;
  };
  
};


PF4PY_NAMESPACE_END

#endif // PF4PY_ELEMSET_H

// Local Variables:
// mode: C++
// End:
