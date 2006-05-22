// $Id: Problem.h,v 1.1.2.8 2006/05/22 18:53:19 dalcinl Exp $

#ifndef PYPF_PROBLEM_H
#define PYPF_PROBLEM_H


#include <string>
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

protected:
  Mesh*   mesh;
  DofMap* dofmap;

public:
  ~Problem();
  Problem(const Problem& problem);
  Problem(Mesh& mesh, DofMap& dofmap);
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

  void buildSolution (Vec state,    Vec solution) const;
  void buildState    (Vec solution, Vec state)    const;

protected:
  virtual void preAssemble(const std::string& jobinfo);
  virtual void baseAssemble(const std::string& jobinfo);
  virtual void postAssemble(const std::string& jobinfo);

};


PYPF_NAMESPACE_END

#endif // PYPF_PROBLEM_H

// Local Variables:
// mode: C++
// End:
