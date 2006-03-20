// -*- c++ -*-
// $Id: Mesh.h,v 1.1.2.4 2006/03/20 16:06:00 rodrigop Exp $

#ifndef PYPF_MESH_H
#define PYPF_MESH_H

#include <string>
#include <vector>
#include "petscfem4py.h"
#include "Object.h"
#include "Nodedata.h"
#include "Elemset.h"

PYPF_NAMESPACE_BEGIN

class Mesh : SMARTPTR(Mesh)
  public Object
{
  friend class Problem;
  
 protected:
  OptionTable* get_opt_table() const; 
  
 protected:
  Nodedata* nodedata;
  std::vector<Elemset*> elemsetlist;

#if !defined(SWIG)
 public:
  Mesh(Mesh::Base*);
#endif

 public:
  ~Mesh();
  Mesh();
  Mesh(const Mesh&);
  Mesh(Nodedata*, const std::vector<Elemset*>&);
  
  Nodedata* getNodedata() const;
  void      setNodedata(Nodedata*);
  
  int      getSize() const;
  Elemset* getElemset(int) const;
  void     setElemset(int, Elemset*);
  void     addElemset(Elemset*);

 public:
  void setUp();
  void clear();
  void view() const;
  

};

PYPF_NAMESPACE_END

#endif // PYPF_MESH_H
