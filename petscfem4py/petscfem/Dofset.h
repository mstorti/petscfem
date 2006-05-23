// $Id: Dofset.h,v 1.1.2.4 2006/05/23 14:51:30 dalcinl Exp $ 

#ifndef PYPF_DOFSET_H
#define PYPF_DOFSET_H

#include <vector>
#include <list>
#include <set>
#include "petscfem4py.h"
#include "Object.h"
#include "Amplitude.h"

PYPF_NAMESPACE_BEGIN

class Dofset :
  public Object
{
  friend class DofMap;

private:
  Dofset();
  
protected:
  struct fixation {
    int node; int field; double value; Amplitude::Base* amp;
    fixation()
      : node(0), field(0), value(0.0), amp(NULL) { }
    fixation(const fixation& F)
      : node(F.node), field(F.field), value(F.value), amp(F.amp) { }
    fixation(int n, int f, double c, Amplitude::Base* a)
      : node(n), field(f), value(c), amp(a) { }
  };
  typedef fixation              Fixation;
  typedef std::vector<Fixation> FixationList;
  FixationList fixations;
  typedef std::set<Amplitude*>  AmplitudeSet;
  AmplitudeSet amplitude;

protected:
  struct constraint {
    int node; int field; double coeff;
    constraint() 
      : node(0), field(0), coeff(0.0) { }
    constraint(const constraint& C)
      : node(C.node), field(C.field), coeff(C.coeff) { }
    constraint(int n, int f, double c)
      : node(n), field(f), coeff(c) { }
  };
  typedef std::vector<constraint> Constraint;
  typedef std::list<Constraint>   ConstraintList;
  ConstraintList constraints;

protected:
  int nnod, ndof;

public:
  ~Dofset();
  Dofset(const Dofset& dofset);
  Dofset(int nnod, int ndof);

 public:
  void addFixations(int n,
		    const int    node[],
		    const int    field[],
		    const double value[]);

  void addFixations(int n,
		    const int    node[],
		    const int    field[],
		    const double value[], 
		    Amplitude*   amplitude);

  void addConstraints(int n,
		      const int node[],
		      const int field[],
		      const double coeff[]);

  void clear();

};

PYPF_NAMESPACE_END

#endif // PYPF_DOFSET_H

// Local Variables:
// mode: C++
// End:
