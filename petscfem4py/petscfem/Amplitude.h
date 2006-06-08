// $Id: Amplitude.h,v 1.1.2.2 2006/06/08 15:44:52 dalcinl Exp $ 

#ifndef PYPF_AMPLITUDE_H
#define PYPF_AMPLITUDE_H

#include "petscfem4py.h"
#include "Object.h"

PYPF_NAMESPACE_BEGIN

class Amplitude : SMARTPTR(Amplitude)
  public Object
{

 private:
  Amplitude(const Amplitude&);

 public:
  enum Kind {
    CONSTANT = 0,
    TEMPORAL = 1,
    NODAL    = 2,
    GENERAL  = 3
  };

 public:
  virtual ~Amplitude();
  Amplitude(Kind kind=GENERAL);

//  public:
//   virtual void init();
//   virtual void clear();

 public:
  virtual double operator()() = 0;
  virtual double operator()(double time) = 0;
  virtual double operator()(int node, int field) = 0;
  virtual double operator()(int node, int field, double time) = 0;

};

PYPF_NAMESPACE_END

#endif // PYPF_AMPLITUDE_H

// Local Variables:
// mode: C++
// End:
