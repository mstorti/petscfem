// $Id$

#include "Amplitude.h"

#include <fem.h>
#include <dofmap.h>

PF4PY_NAMESPACE_BEGIN

typedef Amplitude Amp;

static inline double
mk_time(const TimeData *time_data) {
  if (time_data != NULL)
    return reinterpret_cast<const Time&>(*time_data);
  else
    return 0.0;
}

class AmpProxy : public Amplitude::Impl {
private:
  AmpProxy() : amp(NULL) { }
  AmpProxy(const AmpProxy& A) : amp(A.amp) { }
  AmpProxy& operator=(AmpProxy& A) { this->amp = A.amp; return *this; }
protected:
  Amp*      amp;
public:
  ~AmpProxy() { }
  AmpProxy(Amp* a) : amp(a) { }
  int needs_dof_field_q() 
  { 
    return (this->amp->kind > Amp::TEMPORAL) ? 1 : 0;
  }
  double eval(const TimeData *time_data) 
  {
    switch(this->amp->kind) {
    case Amp::TEMPORAL:
      return this->amp->operator()(mk_time(time_data));
    case Amp::CONSTANT:
      return this->amp->operator()();
    default:
      return 0.0;
    }
  }
  double eval(const TimeData *time_data, int node, int field)
  {
    --node; --field;
    switch(this->amp->kind) {
    case Amp::GENERAL:
      return this->amp->operator()(mk_time(time_data), node, field);
    case Amp::NODAL:
      return this->amp->operator()(node, field);
    default:
      return 0.0;
    }
  }
  void init(::TextHashTable *t) 
  { 
    PF4PY_DELETE_SCLR(t);
  }
  void print() const { }
};

PF4PY_NAMESPACE_END



PF4PY_NAMESPACE_BEGIN

Amplitude::~Amplitude()
{ }

Amplitude::Amplitude(const Amplitude& a)
  : impl(new AmpProxy(this)),
    kind(a.kind)
  { }

Amplitude::Amplitude(Kind kind)
  : impl(new AmpProxy(this)),
    kind(kind)
{ }

Amplitude&
Amplitude::operator=(const Amplitude& A)
{
  this->kind = A.kind;
  return *this;
}

PF4PY_NAMESPACE_END
