// $Id: Problem.h,v 1.1.2.9 2006/06/06 15:45:49 dalcinl Exp $

#ifndef PYPF_PROBLEM_H
#define PYPF_PROBLEM_H


#include <vector>
#include "petscfem4py.h"
#include "Object.h"
#include "Nodeset.h"
#include "Elemset.h"
#include "Mesh.h"
#include "DofMap.h"


PYPF_NAMESPACE_BEGIN

class Problem 
  : public Object
{

private:
  Problem();
  Problem(const Problem& problem);
  Problem(Mesh& mesh, DofMap& dofmap);

protected:
  Mesh*   mesh;
  DofMap* dofmap;
  
public:
  ~Problem();
  Problem(Mesh& mesh,
	  Dofset& dofset);
  Problem(Nodeset& nodeset,
	  const std::vector<Elemset*>& elemsets,
	  Dofset& dofset);
  
  Mesh&   getMesh()   const;
  DofMap& getDofMap() const;

  int  getDim() const;
  int  getSize() const;
  void getSizes(int* nnod, int* ndof) const;

  int  getDofSize() const;
  void getDofSizes(int* local, int* global) const;
  void getDofRange(int* first, int* last)   const;

  void getLocalDofs(int* n, int* dofs[]) const;


  void buildSolution (Vec state, Vec soltn, double t=0.0) const;
  void buildState    (Vec soltn, Vec state) const;

protected:
  void preAssemble();
  void postAssemble();

};


PYPF_NAMESPACE_END

#endif // PYPF_PROBLEM_H

// Local Variables:
// mode: C++
// End:
