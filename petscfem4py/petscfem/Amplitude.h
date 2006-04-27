// $Id: Amplitude.h,v 1.1.2.1 2006/04/27 19:09:17 rodrigop Exp $ 

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
  virtual ~Amplitude();
  Amplitude();

 public:
  virtual double operator()(int node, int field, double time);

};

PYPF_NAMESPACE_END

#endif // PYPF_AMPLITUDE_H

// Local Variables:
// mode: C++
// End:
