// -*- c++ -*-

#ifndef PYPF_MESH_H
#define PYPF_MESH_H

#include <string>
#include <map>
#include "petscfem4py.h"


PYPF_NAMESPACE_BEGIN

PYPF_CLASS(Mesh)
{
  PYPF_CONSTRUCTOR(Mesh)

 public:
  Mesh();
  ~Mesh();
  
  std::string getOption(const std::string& key);
  void setOption(const std::string& name,
		 const std::string& value);

  Nodedata getNodedata();
  void setNodedata(Nodedata& nodedata);

  void addElemset(Elemset& elemset);
  Elemset getElemset(int i);
  int getSize();

  bool hasElemset(const std::string& name);
  Elemset findElemset(const std::string& name);

};


PYPF_NAMESPACE_END

#endif // PYPF_MESH_H
