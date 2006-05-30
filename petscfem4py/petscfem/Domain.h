// $Id: Domain.h,v 1.1.2.4 2006/05/30 20:02:46 dalcinl Exp $

#ifndef PYPF_DOMAIN_H
#define PYPF_DOMAIN_H

#include <vector>
#include "petscfem4py.h"
#include "Object.h"
#include "Nodeset.h"
#include "Elemset.h"
#include "Mesh.h"
#include "DofMap.h"

PYPF_NAMESPACE_BEGIN

class Domain :
  public Object
{

private:
  Domain();

protected:
  Mesh*   mesh;
  Dofset* dofset;
  DofMap* dofmap;
  
public:
  ~Domain();
  Domain(const Domain& domain);
  Domain(Mesh& mesh,
	 Dofset& dofset);
  Domain(Nodeset& nodeset,
	 const std::vector<Elemset*>& elemsets,
	 Dofset& dofset);
  
  // Nodeset&              getNodeset() const;
  // std::vector<Elemset*> getElemset() const;
  // Dofset&               getDofset()  const;
  
  Mesh&   getMesh()   const;
  DofMap& getDofMap() const;

  int  getDim() const;
  int  getSize() const;
  void getSizes(int* nnod, int* ndof) const;

  int  getDofSize() const;
  void getDofSizes(int* local, int* global)  const;
  void getDofRange(int* first, int* last)    const;
  void getDofDist(int* rsize, int* ranges[]) const;

  void getOwnedDofs(int* start, int* end)    const;
  void getGhostDofs(std::vector<int>& gdofs) const;
  void getLocalDofs(std::vector<int>& ldofs) const;

};

PYPF_NAMESPACE_END

#endif // PYPF_DOMAIN_H

// Local Variables:
// mode: C++
// End:
