// $Id: Mesh.h,v 1.1.2.8 2006/04/27 19:09:17 rodrigop Exp $

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
  Nodeset*              nodeset;
  std::vector<Elemset*> elemsetlist;

 public:
  ~Mesh();
  Mesh();
  Mesh(const Mesh& mesh);
  Mesh(Nodeset& nodeset,
       const std::vector<Elemset*>& elemsets);
  
  Nodeset& getNodeset() const;
  void     setNodeset(Nodeset& nodeset);

  Elemset& getElemset(int i) const;
  void     setElemset(int i, Elemset& elemset);
  void     delElemset(int i);
  void     addElemset(Elemset& elemset);
  int      getSize() const;

 public:
  void clear();
  void view() const;

};

PYPF_NAMESPACE_END

#endif // PYPF_MESH_H

// Local Variables:
// mode: C++
// End:
