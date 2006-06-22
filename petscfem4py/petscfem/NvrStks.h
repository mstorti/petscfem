// $Id: NvrStks.h,v 1.1.2.4 2006/06/22 22:35:42 dalcinl Exp $

#ifndef PYPF_NVRSTKS_H
#define PYPF_NVRSTKS_H

#include "petscfem4py.h"
#include "Domain.h"
#include "Application.h"

PYPF_NAMESPACE_BEGIN

class ArgsNS;

class NvrStks : 
  public Application
{
protected:
  NvrStks();

protected:
  double  alpha;
  bool    steady;
  ArgsNS* args;
  
public:
  ~NvrStks();
  NvrStks(const NvrStks& ns);
  NvrStks(Domain& domain, double alpha=1.0, bool steady=false);

  double getAlpha()             { return this->alpha;    };
  bool   getSteady()            { return this->steady;   };
  void   setAlpha(double alpha) { this->alpha  = alpha;  };
  void   setSteady(bool steady) { this->steady = steady; };

  void operator()(double t, Vec x, Vec r, Mat J) const
  { this->assemble(t, x, r, J); }
  void operator()(double t0, Vec x0, double t1, Vec x1, Vec r, Mat J) const
  { this->assemble(t0, x0, t1, x1, r, J); }
  
  void assemble(double t, Vec x, Vec r, Mat J) const;
  void assemble(double t0, Vec x0, double t1, Vec x1, Vec r, Mat J) const;

};

PYPF_NAMESPACE_END

#endif // PYPF_NVRSTKS_H

// Local Variables:
// mode: C++
// End:
