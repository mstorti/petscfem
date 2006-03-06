// -*- c++ -*-
// $Id: Problem.h,v 1.1.2.2 2006/03/06 16:56:04 rodrigop Exp $

#ifndef PYPF_PROBLEM_H
#define PYPF_PROBLEM_H


#include <string>
#include <mpi.h>
#include "petscfem4py.h"
#include "Nodedata.h"
#include "Elemset.h"
#include "Mesh.h"
#include "DofMap.h"


PYPF_NAMESPACE_BEGIN

class Problem
{
protected:

  MPI_Comm  comm;
  int nnod, ndim, ndof;
  Mesh::Base*   mesh;
  DofMap::Base* dofmap;
  bool setupcalled;

public:
  virtual ~Problem();
  Problem();
  Problem(int nnod, int ndim, int ndof);

  void fromFile(const std::string& filename);

  MPI_Comm& getComm();
  Mesh      getMesh();
  DofMap    getDofMap();
  /*
  Nodedata  getNodeData();
  */

  virtual void setUp();

  /*
  virtual void computeResidual() = 0;
  virtual void computeJacobian() = 0;
  */
};


PYPF_NAMESPACE_END

#endif // PYPF_PROBLEM_H
