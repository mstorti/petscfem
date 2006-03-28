// -*- c++ -*-
// $Id: Problem.h,v 1.1.2.4 2006/03/28 22:13:25 rodrigop Exp $

#ifndef PYPF_PROBLEM_H
#define PYPF_PROBLEM_H


#include <string>
#include <vector>
#include <mpi.h>
#include "petscfem4py.h"
#include "Object.h"
#include "Nodedata.h"
#include "Elemset.h"
#include "Mesh.h"
#include "DofMap.h"


PYPF_NAMESPACE_BEGIN

class Problem : 
  public Object
{

private:
  //Problem();

protected:
  int nnod, ndim, ndof;
  Mesh*   mesh;
  DofMap* dofmap;
  bool setupcalled;

public:
  ~Problem();
  Problem();
  Problem(const Problem &);
  Problem(Mesh* mesh, DofMap* dofmap);

  Mesh*   getMesh()   const;
  DofMap* getDofMap() const;

  void read(const std::string& filename);

  void buildState(Vec solution, Vec state);
  void buildSolution(Vec state, Vec solution);

  void setUp();
  void clear();
  

};


PYPF_NAMESPACE_END

#endif // PYPF_PROBLEM_H
