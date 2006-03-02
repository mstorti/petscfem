// -*- c++ -*-
// $Id: Elemset.h,v 1.1.2.3 2006/03/02 21:37:12 rodrigop Exp $

#ifndef PYPF_ELEMSET_H
#define PYPF_ELEMSET_H

#include <string>
#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN

PYPF_CLASS(Elemset) 
{
  PYPF_CTOR(Elemset)

 public:
  Elemset(const std::string& type,
	  const std::string& name="");
  ~Elemset();
  
  void setUp();

  std::string getOption(const std::string& key);
  void        setOption(const std::string& name,
			const std::string& value);

  std::string getType();
  std::string getName();

  void getConnectivity(int* nelem, int* nel, int* icone[]);
  void setConnectivity(int  nelem, int  nel, int  icone[]);

  void getSize(int* nelem, int* nel);

  void setNDof(int ndof);
  int  getNDof();

  void view();

};


PYPF_NAMESPACE_END

#endif // PYPF_ELEMSET_H
