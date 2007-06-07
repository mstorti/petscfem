// $Id$ 

#ifndef PF4PY_AMPLITUDE_H
#define PF4PY_AMPLITUDE_H

#include <memory>
#include "petscfem4py.h"

#include "Object.h"

PF4PY_NAMESPACE_BEGIN

class Amplitude
  : public Object
{
  friend class Dofset;
  friend class Domain;

public:
  typedef ::Amplitude Impl;

public:
  enum Kind {
    CONSTANT = 0,
    TEMPORAL = 1,
    NODAL    = 2,
    GENERAL  = 3
  };

#if !defined(SWIG)
protected:
  std::auto_ptr< ::Amplitude > impl;
public: 
  inline operator ::Amplitude*() const { return this->impl.get(); }
#endif

public:
  Kind kind;

public:
  virtual ~Amplitude();
protected:
  Amplitude(const Amplitude&);
  Amplitude& operator=(const Amplitude&);
public:
  Amplitude(Kind kind=GENERAL);
  
public:
  virtual double operator()() = 0;
  virtual double operator()(double time) = 0;
  virtual double operator()(int node, int field) = 0;
  virtual double operator()(double time, int node, int field) = 0;

};

PF4PY_NAMESPACE_END

#endif // PF4PY_AMPLITUDE_H

// Local Variables:
// mode: C++
// End:
