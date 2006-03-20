// -*- c++ -*-
// $Id: Problem.h,v 1.1.2.3 2006/03/20 16:06:00 rodrigop Exp $

#ifndef PYPF_PROBLEM_H
#define PYPF_PROBLEM_H


#include <string>
#include <map>
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

protected:
  virtual OptionTable* get_opt_table() const;

private:
  //Problem();
  Problem(const Problem &);

protected:
  MPI_Comm comm;
  int nnod, ndim, ndof;
  Mesh*   mesh;
  DofMap* dofmap;
  bool setupcalled;

public:
  ~Problem();
  Problem(int nnod, int ndim, int ndof);

  Mesh*      getMesh()   const;
  DofMap*    getDofMap() const;

  Problem();
  void fromFile(const std::string& filename);

  void buildState(Vec solution, Vec state);
  void buildSolution(Vec state, Vec solution);

  void setUp();
  void clear();
  

};


PYPF_NAMESPACE_END

#endif // PYPF_PROBLEM_H
