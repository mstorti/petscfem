// -*- c++ -*-
// $Id: Problem.h,v 1.1.2.6 2006/03/30 15:40:05 rodrigop Exp $

#ifndef PYPF_PROBLEM_H
#define PYPF_PROBLEM_H


#include <string>
#include <vector>
#include <mpi.h>
#include "petscfem4py.h"
#include "Object.h"
#include "Nodeset.h"
#include "Elemset.h"
#include "Mesh.h"
#include "DofMap.h"


PYPF_NAMESPACE_BEGIN

class Problem : 
  public Object
{

private:
  Problem();
  Problem(const Problem&);

protected:
  int nnod, ndim, ndof;
  Mesh*   mesh;
  DofMap* dofmap;
  bool setupcalled;

public:
  ~Problem();
  Problem(Mesh*, DofMap*);

  void setUp();
  
  Mesh*   getMesh()   const;
  DofMap* getDofMap() const;

  void getDofSizes (int* local, int* global) const;
  void getDofRange (int* first, int* last)   const;

  void buildSolution (Vec state,    Vec solution) const;
  void buildState    (Vec solution, Vec state)    const;

  //void read(const std::string& filename);
};


PYPF_NAMESPACE_END

#endif // PYPF_PROBLEM_H
