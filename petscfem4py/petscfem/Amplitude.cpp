// $Id: Amplitude.cpp,v 1.1.2.2 2006/06/05 20:38:57 dalcinl Exp $

#include "Amplitude.h"

#include <fem.h>
#include <dofmap.h>

class AmpCaller : public Amplitude {
private:
  AmpCaller() : amp(NULL) { }
protected:
  PYPF_NAMESPACE::Amplitude* amp;
public:
  ~AmpCaller() { }
  AmpCaller(PYPF_NAMESPACE::Amplitude* amp) : amp(amp) { }
  double eval(const TimeData *time_data, int node, int field) {
    double time = 0.0;
    if (time_data) time = reinterpret_cast<const Time&>(*time_data);
    return this->amp->operator()(node-1, field-1, time);
  }
public:
  double eval(const TimeData *time_data) { return 1.0; }
  void   init(::TextHashTable *t) { PYPF_DELETE_SCLR(t); }
  void   print() const { }
  int    needs_dof_field_q() { return 1; }
};


PYPF_NAMESPACE_BEGIN

Amplitude::~Amplitude()
{
  Amplitude::Base* amp = *this;
  if (amp == NULL) return;
  delete amp;
}

Amplitude::Amplitude()
  : Handle(new AmpCaller(this)), 
    Object()
{ }

Amplitude::Amplitude(const Amplitude& amp)
  : Handle(new AmpCaller(this)),
    Object(amp)
{ }

double
Amplitude::operator()(int node, int field, double time)
{
  return 1.0;
}


PYPF_NAMESPACE_END
