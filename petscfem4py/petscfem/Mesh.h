// -*- c++ -*-
// $Id: Mesh.h,v 1.1.2.3 2006/03/06 16:56:04 rodrigop Exp $

#ifndef PYPF_MESH_H
#define PYPF_MESH_H

#include <string>
#include <vector>
#include "petscfem4py.h"
#include "Nodedata.h"
#include "Elemset.h"

PYPF_NAMESPACE_BEGIN

PYPF_CLASS(Mesh)
{
  PYPF_CTOR_FROM_PTR(Mesh)
  PYPF_OBJ_GETOPTTBL_DECL

 protected:

 public:
  ~Mesh();
  Mesh();
  
//   std::string getOption(const std::string& key);
//   void        setOption(const std::string& key,
// 			const std::string& value);

  Nodedata getNodeData();
  
  Elemset  getElemset(int i);
  Elemset  getElemset(const std::string& name);
  
  Elemset  addElemset(const std::string& name,
		      const std::string& type);

//   void    addElemset(Elemset& elemset);
//   bool    hasElemset(const std::string& name);
//   Elemset findElemset(const std::string& name);

  int getSize();

  friend class Problem;
};

PYPF_NAMESPACE_END

#endif // PYPF_MESH_H
