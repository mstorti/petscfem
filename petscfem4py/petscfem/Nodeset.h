// $Id: Nodeset.h,v 1.1.2.9 2006/06/07 16:27:16 dalcinl Exp $

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
  Nodeset();
  void chk_sizes(int nnod, int ndim, int nval) const;
  void set_sizes(int nnod, int ndim, int nval);
  void set_array(const double array[]);
  void touch() const;

 protected:
  int nnod, ndim, nval;
  std::vector<double> nodedata;

 public:
  ~Nodeset();
  Nodeset(const Nodeset& nodeset);
  Nodeset(int nnod, int ndim, int nval=0);
  Nodeset(int nnod, int ndim, const double xnod[]);

  void getSizes(int* nnod, int* ndim, int* nval) const;
  void setSizes(int nnod, int ndim, int nval=0);
  void getArray(const double* array[]) const;
  void setArray(const double  array[]);

  int  getDim()  const;
  int  getSize() const;
  void getNode(int i, int* n, const double* node[]) const;
  void setNode(int i, int  n, const double  node[]);

 public:
  void clear();
  void view() const;

};

PYPF_NAMESPACE_END

#endif // PYPF_NODESET_H

// Local Variables:
// mode: C++
// End:
