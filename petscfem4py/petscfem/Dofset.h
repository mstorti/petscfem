// $Id: Dofset.h,v 1.1.2.8 2006/06/22 22:34:59 dalcinl Exp $ 

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

protected:
  struct AmpSet : private std::set<Amplitude*> {
    using std::set<Amplitude*>::iterator;
    using std::set<Amplitude*>::const_iterator;
    using std::set<Amplitude*>::begin;
    using std::set<Amplitude*>::end;
    ~AmpSet()  { this->clear(); }
    AmpSet() : std::set<Amplitude*>() { }
    AmpSet(const AmpSet& as) : std::set<Amplitude*>(as) {
      iterator s = this->begin();
      while (s != this->end()) { Amplitude* a = *s++; PYPF_INCREF(a); }
    }
    void add(Amplitude* a) {
      if (a == NULL) return;
      std::pair<iterator,bool> p = std::set<Amplitude*>::insert(a);
      if (p.second) PYPF_INCREF(a);
    }
    void clear() {
      iterator s = this->begin();
      while (s != this->end()) { Amplitude* a = *s++; PYPF_DECREF(a); }
      std::set<Amplitude*>::clear();
    }
  };
  typedef AmpSet AmplitudeSet;
  AmplitudeSet amplitudes;

protected:
  struct constraint {
    double coeff; int node; int field;
    constraint() 
      : coeff(0.0), node(0), field(0) { }
    constraint(const constraint& C)
      : coeff(C.coeff), node(C.node), field(C.field) { }
    constraint(double c, int n, int f)
      : coeff(c), node(n), field(f) { }
  };
  typedef std::vector<constraint> Constraint;
  typedef std::list<Constraint>   ConstraintList;
  ConstraintList constraints;


protected:
  Dofset();

protected:
  inline void chk_fixa(int, int);
  inline void chk_fixa(int, const int[], const int[]);

  inline void chk_fixa(int, int, double);
  inline void chk_fixa(int, const int[], const int[], const double[]);

  inline void add_fixa(int, int, double, 
		       Amplitude::Base* a=NULL);
  inline void add_fixa(int, const int[], const int[], const double[],
		       Amplitude::Base* a=NULL);
 
protected:
  int nnod, ndof;

public:
  ~Dofset();
  Dofset(const Dofset& dofset);
  Dofset(int nnod, int ndof);
  Dofset(int nnod, int ndof, MPI_Comm comm);

 public:

  void getSizes(int* nnod, int* ndof) const;

  void addFixations(int n,
		    const int    node[],
		    const int    field[],
		    const double value[],
		    Amplitude* amplitude=NULL);

  void addConstraints(int n,
		      const double coeff[],
		      const int    node[],
		      const int    field[]);
public:
  void view() const;
  void clear();

};

PYPF_NAMESPACE_END

#endif // PYPF_DOFSET_H

// Local Variables:
// mode: C++
// End:
