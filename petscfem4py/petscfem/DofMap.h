// $Id: DofMap.h,v 1.1.2.9 2006/06/08 15:44:52 dalcinl Exp $

#ifndef PYPF_DOFMAP_H
#define PYPF_DOFMAP_H

#include <set>
#include <mpi.h>
#include "petscfem4py.h"
#include "Object.h"
#include "Mesh.h"
#include "Amplitude.h"
#include "Dofset.h"

PYPF_NAMESPACE_BEGIN


class DofMap : SMARTPTR(DofMap)
  public Object
{
  friend class Problem;

 private:
  DofMap();
  DofMap(const DofMap& dofmap);

 protected:
  Dofset::AmplitudeSet ampset;

 protected:
  void build(const Mesh& mesh, const Dofset& dofset);
  void add_fixation  (const Dofset::Fixation&   fixation);
  void add_constraint(const Dofset::Constraint& constraint);
  
 public:
  ~DofMap();
  DofMap(Mesh& mesh, Dofset& dofset);
  DofMap(Mesh& mesh, Dofset& dofset, MPI_Comm comm);

  int  getSize() const;
  int  getLocalSize() const;
  void getSizes(int* local, int* global) const;
  void getRange(int* first, int* last) const;
  void getDist(int* rsize, int* ranges[]) const;

  int getNNod() const;
  int getNDof() const;
  
  void apply(int nstt, const double stt[],
	     int nsol, double sol[],
	     double time=0.0) const;
  void solve(int nsol, const double sol[],
	     int nstt, double stt[]) const;

 public:
  void view() const;

};


PYPF_NAMESPACE_END

#endif // PYPF_DOFMAP_H

// Local Variables:
// mode: C++
// End:
