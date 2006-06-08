// $Id: Amplitude.cpp,v 1.1.2.3 2006/06/08 15:44:52 dalcinl Exp $

#include "Amplitude.h"

#include <fem.h>
#include <dofmap.h>

typedef PYPF_NAMESPACE::Amplitude Amp;

PYPF_NAMESPACE_BEGIN

static inline double get_time(const TimeData *time_data) {
  if (time_data != NULL)
    return reinterpret_cast<const Time&>(*time_data);
  else
    return 0.0;
}

class AmpProxy : public Amp::Base {
protected:
  AmpProxy() : amp(NULL), kind(Amp::GENERAL) { }
  Amp*      amp;
  Amp::Kind kind;
public:
  AmpProxy(const AmpProxy& A) : amp(A.amp), kind(A.kind) { }
  AmpProxy(Amp* a, Amp::Kind k) : amp(a), kind(k) { }
  int needs_dof_field_q() { 
    return this->kind > Amp::TEMPORAL ? 1 : 0;
  }
  double eval(const TimeData *time_data) {
    switch(this->kind) {
    case Amp::TEMPORAL:
      return this->amp->operator()(get_time(time_data));
    case Amp::CONSTANT:
      return this->amp->operator()();
    default:
      return 1.0;
    }
  }
  double eval(const TimeData *time_data, int node, int field) {
    --node; --field;
    switch(this->kind) {
    case Amp::GENERAL:
      return this->amp->operator()(node, field, get_time(time_data));
    case Amp::NODAL:
      return this->amp->operator()(node, field);
    default:
      return 1.0;
    }
  }
  void  init(::TextHashTable *t) { PYPF_DELETE_SCLR(t); }
  void  print() const { }
};

PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN

Amplitude::~Amplitude()
{
  Amplitude::Base* amp = *this;
  if (amp) delete amp;
}

Amplitude::Amplitude(Kind kind)
  : Handle(new AmpProxy(this, kind)),
    Object()
{ }

PYPF_NAMESPACE_END
