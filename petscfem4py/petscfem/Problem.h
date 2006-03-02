// -*- c++ -*-
// $Id: Problem.h,v 1.1.2.1 2006/03/02 21:37:12 rodrigop Exp $

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

  int nnod;
  int ndim;
  int ndof;

  MPI_Comm  comm;
  ::Mesh*   mesh;
  ::DofMap* dofmap;

  bool setupcalled;

public:
  Problem();
  virtual ~Problem();
  
  Problem(int nnod, int ndim, int ndof);

  void fromFile(const std::string& filename);

  MPI_Comm* getComm();
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
