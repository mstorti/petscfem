// -*- c++ -*-
// $Id: DofMap.h,v 1.1.2.6 2006/03/30 15:18:14 rodrigop Exp $

#ifndef PYPF_DOFMAP_H
#define PYPF_DOFMAP_H

#include <mpi.h>
#include "petscfem4py.h"
#include "Object.h"
#include "Mesh.h"

PYPF_NAMESPACE_BEGIN


class DofMap : SMARTPTR(DofMap)
  public Object
{
  friend class Problem;

#if !defined(SWIG)
 public:
  DofMap(DofMap::Base*);
#endif

 private:
  DofMap();
  DofMap(const DofMap&);
  
 protected:
  int  nnod, ndof;
  bool frozen;

 public:
  ~DofMap();
  DofMap(Mesh*, int ndof);

  void addFixations   (int n, int node[], int field[], double value[]);
  void addConstraints (int n, int node[], int field[], double coeff[]); 

  void getSizes  (int* local, int* global)   const;
  void getRange  (int* start, int* end)      const;
  void getRanges (int* rsize, int* ranges[]) const;

  
  int getNNod() const;
  int getNDof() const;
  int getNFix() const;
  
  void setUp();
  void view() const;

};


PYPF_NAMESPACE_END

#endif // PYPF_DOFMAP_H
