// -*- c++ -*-
// $Id: Mesh.h,v 1.1.2.6 2006/03/30 15:18:14 rodrigop Exp $

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
  friend class DofMap;
  friend class Problem;
  
#if !defined(SWIG)
 public:
  Mesh(Mesh::Base*);
#endif

 protected:
  Nodedata*             nodedata;
  std::vector<Elemset*> elemsetlist;

 public:
  ~Mesh();
  Mesh();
  Mesh(const Mesh&);
  Mesh(Nodedata*, const std::vector<Elemset*>&);
  
  Nodedata* getNodedata() const;
  void      setNodedata(Nodedata*);

  int       getSize() const;
  Elemset*  getElemset(int) const;
  void      setElemset(int, Elemset*);
  void      addElemset(Elemset*);

 public:
  void setUp();
  void clear();
  void view() const;

};

PYPF_NAMESPACE_END

#endif // PYPF_MESH_H
