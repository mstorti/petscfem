// $Id: Nodeset.h,v 1.1.2.8 2006/06/06 15:15:33 dalcinl Exp $

#ifndef PYPF_NODESET_H
#define PYPF_NODESET_H

#include <vector>
#include "petscfem4py.h"
#include "Object.h"

PYPF_NAMESPACE_BEGIN

class Nodeset : SMARTPTR(Nodeset)
  public Object
{

 private:
  Nodeset();
  
 protected:
  int nnod, ndim, nval;
  std::vector<double> nodedata;

 protected:
  void chk_sizes(int nnod, int ndim, int nval) const;
  void set_sizes(int nnod, int ndim, int nval);
  void set_array(const double array[]);
  void touch() const;

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
  void sync(int root = 0);
  void clear();
  void view() const;

};

PYPF_NAMESPACE_END

#endif // PYPF_NODESET_H

// Local Variables:
// mode: C++
// End:
