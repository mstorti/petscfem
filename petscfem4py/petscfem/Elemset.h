// -*- c++ -*-
// $Id: Elemset.h,v 1.1.2.7 2006/03/30 15:18:14 rodrigop Exp $

#ifndef PYPF_ELEMSET_H
#define PYPF_ELEMSET_H

#include <string>
#include <vector>
#include "petscfem4py.h"
#include "Object.h"


PYPF_NAMESPACE_BEGIN

class Elemset : SMARTPTR(Elemset)  
  public Object
{
  friend class Mesh;
  friend class DofMap;
  friend class Problem;

#if !defined(SWIG)
 public:
  Elemset(Elemset::Base*);
#endif
  
 private:
  Elemset();
  
 protected:
  int nelem, nel;
  std::vector<int> icone;

 public:
  ~Elemset();
  Elemset(const Elemset&);
  Elemset(const std::string& type, const std::string& name="");
  
  std::string getType() const;
  std::string getName() const;

  void getSize(int* nelem, int* nel) const;
  void getConnectivity(int* nelem, int* nel, int* icone[]) const;
  void setConnectivity(int  nelem, int  nel, int  icone[]);

  typedef std::vector<int> Elem;
  Elem getElem(int i) const;
  void setElem(int i, const Elem& elem);

  int  getNDof() const;
  void setNDof(int ndof);

  void setUp();
  void clear();
  void view() const;

};

PYPF_NAMESPACE_END

#endif // PYPF_ELEMSET_H
