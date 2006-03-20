// -*- c++ -*-
// $Id: Elemset.h,v 1.1.2.5 2006/03/20 16:06:00 rodrigop Exp $

#ifndef PYPF_ELEMSET_H
#define PYPF_ELEMSET_H

#include <string>
#include "petscfem4py.h"
#include "Object.h"


PYPF_NAMESPACE_BEGIN

class Elemset : SMARTPTR(Elemset)
  public Object
{
  friend class Mesh;

protected: 
  OptionTable* get_opt_table() const;
  
private:
  Elemset();
  
#if !defined(SWIG)
public:
  Elemset(Elemset::Base*);
#endif
  
public:
  ~Elemset();
  Elemset(const Elemset&);
  Elemset(const std::string& type,
	  const std::string& name="");
  
  std::string getType() const;
  std::string getName() const;

  void getSize(int* nelem, int* nel) const;
  void getConnectivity(int* nelem, int* nel, int* icone[]) const;
  void setConnectivity(int  nelem, int  nel, int  icone[]);

  int  getNDof() const;
  void setNDof(int ndof);

  void setUp();
  void clear();
  void view() const;

};

PYPF_NAMESPACE_END

#endif // PYPF_ELEMSET_H
