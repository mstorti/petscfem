// -*- c++ -*-

#ifndef PYPF_NODEDATA_H
#define PYPF_NODEDATA_H

#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN

PYPF_CLASS(Nodedata)
{

  PYPF_CONSTRUCTOR(Nodedata)
    
 public:
  
  Nodedata(int nnod, int ndim, double xyz[]);
  ~Nodedata();
  /*
    Nodedata(int nnod, int ndim, const double xyz[],
    int nfields, const double fielddata[]);
  */

  std::string getOption(const std::string& key);
  void setOption(const std::string& name,
		 const std::string& value);

  int  getSize();
  int  getDim();
  void getArray(int* nnod, int* ndim, double* xyz[]);
  void setArray(int  nnod, int  ndim, double  xyz[]);

  /*  
  view();
  */
};


PYPF_NAMESPACE_END

#endif // PYPF_NODEDATA_H
