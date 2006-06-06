// $Id: Mesh.h,v 1.1.2.9 2006/06/06 15:44:27 dalcinl Exp $

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
 private:
  Mesh();
  
 protected:
  Nodeset*              nodeset;
  std::vector<Elemset*> elemsetlist;

 public:
  ~Mesh();
  Mesh(const Mesh& mesh);
  Mesh(Nodeset& nodeset, const std::vector<Elemset*>& elemsets);
  
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
