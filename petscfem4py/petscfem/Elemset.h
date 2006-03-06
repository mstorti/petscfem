// -*- c++ -*-
// $Id: Elemset.h,v 1.1.2.4 2006/03/06 16:56:04 rodrigop Exp $

#ifndef PYPF_ELEMSET_H
#define PYPF_ELEMSET_H

#include <string>
#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN

PYPF_CLASS(Elemset) 
{
  PYPF_CTOR_FROM_PTR(Elemset)
  PYPF_OBJ_GETOPTTBL_DECL

 public:
  ~Elemset();
  Elemset();

  Elemset(const std::string& type,
	  const std::string& name="");
  
  void setUp();

  std::string getType();
  std::string getName();

  void getConnectivity(int* nelem, int* nel, int* icone[]);
  void setConnectivity(int  nelem, int  nel, int  icone[]);

  void getSize(int* nelem, int* nel);

  void setNDof(int ndof);
  int  getNDof();

  void view();

  friend class Mesh;
};


PYPF_NAMESPACE_END

#endif // PYPF_ELEMSET_H
