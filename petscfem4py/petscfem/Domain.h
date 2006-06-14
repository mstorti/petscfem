// $Id: Domain.h,v 1.1.2.7 2006/06/14 19:10:08 dalcinl Exp $

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
  Dofset* dofset;
  Mesh*   mesh;
  DofMap* dofmap;
  
public:
  ~Domain();
  Domain(const Domain& domain);
  Domain(Nodeset& nodeset,
	 const std::vector<Elemset*>& elemsets,
	 Dofset& dofset);
  Domain(Nodeset& nodeset,
	 const std::vector<Elemset*>& elemsets,
	 Dofset& dofset, MPI_Comm comm);
  
  // std::vector<Elemset*> getElemset() const;
  
  Nodeset& getNodeset() const;
  Dofset&  getDofset()  const;
  Mesh&    getMesh()    const;
  DofMap&  getDofMap()  const;

  int  getDim() const;
  void getSizes(int* nnod, int* ndof) const;

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
