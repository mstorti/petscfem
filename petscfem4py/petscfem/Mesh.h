// -*- c++ -*-
// $Id: Mesh.h,v 1.1.2.2 2006/03/02 21:37:12 rodrigop Exp $

#ifndef PYPF_MESH_H
#define PYPF_MESH_H

#include <string>
#include "petscfem4py.h"
#include "Nodedata.h"
#include "Elemset.h"

PYPF_NAMESPACE_BEGIN

PYPF_CLASS(Mesh)
{
  PYPF_CTOR(Mesh)

 public:
  Mesh();
  ~Mesh();
  
  std::string getOption(const std::string& key);
  void        setOption(const std::string& key,
			const std::string& value);

  Nodedata getNodeData();
  
  Elemset  getElemset(int i);
  Elemset  getElemset(const std::string& name);
  
  Elemset  addElemset(const std::string& name,
		      const std::string& type);

//   void    addElemset(Elemset& elemset);
//   bool    hasElemset(const std::string& name);
//   Elemset findElemset(const std::string& name);

  int getSize();

};

PYPF_NAMESPACE_END

#endif // PYPF_MESH_H
