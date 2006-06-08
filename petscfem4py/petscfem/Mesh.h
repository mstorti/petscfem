// $Id: Mesh.h,v 1.1.2.10 2006/06/08 15:44:52 dalcinl Exp $

#ifndef PYPF_MESH_H
#define PYPF_MESH_H

#include <vector>
#include "petscfem4py.h"
#include "Object.h"
#include "Nodeset.h"
#include "Elemset.h"

PYPF_NAMESPACE_BEGIN

class Mesh : SMARTPTR(Mesh)
  public Object
{
 protected:
  Mesh();
  
 protected:
  Nodeset*              nodeset;
  std::vector<Elemset*> elemsetlist;

 public:
  ~Mesh();
  Mesh(const Mesh& mesh);
  Mesh(Nodeset& nodeset,
       const std::vector<Elemset*>& elemsets);
  Mesh(Nodeset& nodeset,
       const std::vector<Elemset*>& elemsets, 
       MPI_Comm comm);
  
  Nodeset& getNodeset() const;
  Elemset& getElemset(int i) const;
  int      getSize() const;

 public:
  void view() const;

};

PYPF_NAMESPACE_END

#endif // PYPF_MESH_H

// Local Variables:
// mode: C++
// End:
