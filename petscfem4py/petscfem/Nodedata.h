// -*- c++ -*-

#ifndef PYPF_NODEDATA_H
#define PYPF_NODEDATA_H

#include <string>
#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN

PYPF_CLASS(Nodedata)
{

  PYPF_CTOR(Nodedata)
    
 public:
  
  Nodedata();
  ~Nodedata();

  std::string getOption(const std::string& key);
  void        setOption(const std::string& key,
			const std::string& value);

  void getData(int* nnod, int* ndim, double* xyz[]);
  void setData(int  nnod, int  ndim, double  xyz[]);
  void getSize(int* nnod, int* ndim);

  void view();
};


PYPF_NAMESPACE_END

#endif // PYPF_NODEDATA_H
