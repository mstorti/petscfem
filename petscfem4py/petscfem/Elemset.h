// -*- c++ -*-

#ifndef PYPF_ELEMSET_H
#define PYPF_ELEMSET_H

#include <string>
#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN

PYPF_CLASS(Elemset) 
{
  PYPF_CONSTRUCTOR(Elemset)

 public:
  Elemset(const std::string& type,
	  const std::string& name="");
  ~Elemset();
  
  std::string getOption(const std::string& key);
  void setOption(const std::string& name,
		 const std::string& value);
  void setUp();

  std::string getType();
  std::string getName();

  void getConnectivity(int* nelem, int* nel, int* icone[]);
  void setConnectivity(int  nelem, int  nel, int  icone[]);
  int  getSize();
  
  void setNDof(int ndof);
  int  getNDof();

};


PYPF_NAMESPACE_END

#endif // PYPF_NODEDATA_H
