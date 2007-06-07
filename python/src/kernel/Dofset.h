// $Id$ 

#ifndef PF4PY_DOFSET_H
#define PF4PY_DOFSET_H

#include <memory>
#include <vector>
#include <list>
#include <set>

#include "petscfem4py.h"

#include "Comm.h"
#include "Object.h"
#include "Amplitude.h"

PF4PY_NAMESPACE_BEGIN

class Dofset
  : public Object
{
  friend class Mesh;
  friend class Domain;

public:
  typedef ::Dofmap Impl;

protected:
  void add_fixa(int, int, double, Amplitude* a=0);
  void add_fixations(int n, 
		     const int[],const int[], const double[],
		     Amplitude* a=0);
  void add_constraints(int n, const double[],const int[],const int[]);
  void del_fixations();
  void del_constraints();
  void clear();

protected:
  struct fixation {
    int node; int field; double value; Amplitude* amp;
    fixation()
      : node(-1), field(-1), value(0.0), amp(NULL) { }
    fixation(const fixation& F)
      : node(F.node), field(F.field), value(F.value), amp(F.amp) { }
    fixation(int n, int f, double c, Amplitude* a=0)
      : node(n), field(f), value(c), amp(a) { }
  };
  typedef fixation              Fixation;
  typedef std::vector<Fixation> FixationList;

protected:
  struct constraint {
    double coeff; int node; int field;
    constraint()
      : coeff(0.0), node(-1), field(-1) { }
    constraint(const constraint& C)
      : coeff(C.coeff), node(C.node), field(C.field) { }
    constraint(double c, int n, int f)
      : coeff(c), node(n), field(f) { }
  };
  typedef std::vector<constraint> Constraint;
  typedef std::list<Constraint>   ConstraintList;

protected:
  MPI_Comm           comm;

  int                nnod;
  int                ndof;
  
  std::vector<float> tpwgts;
  std::vector<int>   ownership;
  std::vector<int>   ghosts;

  FixationList       fixations;
  RefSet<Amplitude>  amplitudes;
  ConstraintList     constraints;

#if !defined(SWIG)
protected:
  class Proxy; friend class Proxy;
  std::auto_ptr<Proxy> proxy;
  Dofset::Impl* getimpl() const;
public:
  inline operator Dofset::Impl*() const { return this->getimpl(); }
#endif

protected:
  void setup(Domain*);

public:
  ~Dofset();
private:
  Dofset();
  Dofset(const Dofset&);
  Dofset& operator=(const Dofset&);
protected:
  Dofset(MPI_Comm comm, int nnod, int ndof);

public:
  Comm getComm() const { return this->comm; }
  int  getNNod() const { return this->nnod; }
  int  getNDof() const { return this->ndof; }

public:
  const std::vector<float>& getWeights() const;
  void setWeights(const std::vector<float>& weights);
  
  std::pair<int,int> getSizes() const;
  std::pair<int,int> getSizes(int rank) const;
  std::pair<int,int> getRange() const;
  std::pair<int,int> getRange(int rank) const;
  std::vector<int>   getDist() const;

public:
  void getGhostDofs(std::vector<int>& gdofs) const;
  void getLocalDofs(std::vector<int>& ldofs) const;
  void getFieldDofs(int field, std::vector<int>& fdofs) const;

};

PF4PY_NAMESPACE_END

#endif // PF4PY_DOFSET_H

// Local Variables:
// mode: C++
// End:
