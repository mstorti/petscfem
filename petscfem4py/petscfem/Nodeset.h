// $Id: Nodeset.h,v 1.1.2.7 2006/06/05 23:54:01 dalcinl Exp $

#ifndef PYPF_NODESET_H
#define PYPF_NODESET_H

#include <vector>
#include "petscfem4py.h"
#include "Object.h"

PYPF_NAMESPACE_BEGIN

class Nodeset : SMARTPTR(Nodeset)
  public Object
{

 protected:
  int nnod, ndim, nval;
  std::vector<double> nodedata;

 public:
  ~Nodeset();
  Nodeset();
  Nodeset(const Nodeset& nodeset);
  Nodeset(int nnod, int ndim, const double xnod[]);
  Nodeset(int ndim, int nnod, int nval, const double data[]);

  int  getDim() const;
  void setDim(int ndim);

  void getDataSize(int* nnod, int* nval) const;
  void getData(int* nnod, int* nval, const double* data[]) const;
  void setData(int  nnod, int  nval, const double  data[]);

  int  getSize() const;
  void getNode(int i, int* n, const double* node[]) const;
  void setNode(int i, int  n, const double  node[]);

 public:
  void sync(int root = 0);
  void clear();
  void view() const;

};

PYPF_NAMESPACE_END

#endif // PYPF_NODESET_H

// Local Variables:
// mode: C++
// End:
