// $Id: Nodeset.h,v 1.1.2.4 2006/06/05 22:14:02 dalcinl Exp $

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
  int ndim;
  int nnod, nval;
  std::vector<double> nodedata;

 public:
  ~Nodeset();
  Nodeset();
  Nodeset(const Nodeset& nodeset);
  Nodeset(int nnod, int nval, const double data[]);
  Nodeset(int ndim, int nnod, int nval, const double data[]);

  int  getDim() const;
  void setDim(int ndim);

  void getData(int* nnod, int* nval, const double* data[]) const;
  void setData(int  nnod, int  nval, const double  data[]);
  void getDataSize(int* nnod, int* nval) const;

  void getNode(int i, int* n, const double* data[]) const;
  void setNode(int i, int  n, const double  data[]);
  int  getSize() const;

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
