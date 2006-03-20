// -*- c++ -*-
// $Id: DofMap.h,v 1.1.2.4 2006/03/20 16:06:00 rodrigop Exp $

#ifndef PYPF_DOFMAP_H
#define PYPF_DOFMAP_H

#include <mpi.h>
#include "petscfem4py.h"
#include "Object.h"

PYPF_NAMESPACE_BEGIN


class DofMap : SMARTPTR(DofMap)
  public Object
{

protected:
  virtual OptionTable* get_opt_table() const;

 private:
  DofMap(const DofMap&);
  
#if !defined(SWIG)
 public:
  DofMap(DofMap::Base*);
#endif

 public:
  DofMap();
  ~DofMap();
  DofMap(int nnod, int ndof);

  void addFixations(int n, int node[], int field[], double value[]);
  void addConstraints(int n, int node[], int field[], double coeff[]);

  void  getSizes(int* local, int* global)    const;
  void  getRange(int* start, int* end)       const;
  void  getRanges(int* rsize, int* ranges[]) const;

  
  void getNnod(int* nnod);
  void getFixSize(int* neq_fix);
  
  void setUp();
  void view() const;

};


PYPF_NAMESPACE_END

#endif // PYPF_DOFMAP_H
