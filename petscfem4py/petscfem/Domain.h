// $Id: Domain.h,v 1.1.2.1 2006/05/24 21:04:02 dalcinl Exp $

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
  Domain(const Domain& domain);

protected:
  Mesh*   mesh;
  Dofset* dofset;
  DofMap* dofmap;
  
public:
  ~Domain();
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
  void getDofSizes(int* local, int* global) const;
  void getDofRange(int* first, int* last)   const;

  void getLocalDofs(std::vector<int>& ldofs) const;
  void getDofGraph(std::vector<int>& xadj,
		   std::vector<int>& adjncy) const;

};

PYPF_NAMESPACE_END

#endif // PYPF_DOMAIN_H

// Local Variables:
// mode: C++
// End:
