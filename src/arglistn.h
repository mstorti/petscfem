// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: arglistn.h,v 1.3 2001/05/30 03:58:50 mstorti Exp $

#ifndef ARGLISTN_H
#define ARGLISTN_H

#include <vector>
#include <string>
#include <mat.h>
#include <vec.h>
#include "libretto.h"

class ArgEntry {
public:
  virtual void chunk_pre_operations(int chunk_size,int ndoft) {};
};

class OutArg : public ArgEntry {
public:
  virtual void chunk_pre_operations(int chunk_size,int ndoft);
private:
  double *retval;
}

class OutVectorArg : public ArgEntry {
public:
  OutVectorArg(Vec &x);
  ~OutVectorArg();
  void chunk_pre_operations(int chunk_size,int ndoft);
private:
  double *retval;
}

typedef vector<ArgEntry *> ArgList;

#endif
