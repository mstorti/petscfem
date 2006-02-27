#include "Nodedata.h"

#include <fem.h>

//PyPF::Nodedata::Nodedata() : Ptr() { }
//PyPF::Nodedata::Nodedata(const PyPF::Nodedata & nd) : Ptr(nd) { }

PyPF::Nodedata::Nodedata(int nnod, int ndim, double xyz[])
  : Ptr(new ::Nodedata)
{ 
  (*this)->nnod = nnod;
  (*this)->ndim = ndim;
  (*this)->nu   = ndim;
  (*this)->nodedata = new double[nnod*ndim];
  this->setArray(nnod, ndim, xyz);

  (*this)->options = new TextHashTable;
}

PyPF::Nodedata::~Nodedata()
{ 
  /* delete[] this->ptr; */
}

std::string
PyPF::Nodedata::getOption(const std::string& key)
{
  const char* value = NULL;
  (*this)->options->get_entry(key.c_str(), value);
  return value;
}

void 
PyPF::Nodedata::setOption(const std::string& key,
			  const std::string& value)
{
  (*this)->options->add_entry(key.c_str(), value.c_str());
}


int
PyPF::Nodedata::getSize()
{
  return (*this)->nnod;
}

int
PyPF::Nodedata::getDim() 
{
  return (*this)->ndim;
}

void
PyPF::Nodedata::getArray(int* nnod, int* ndim, double* xyz[])
{
  *nnod = (*this)->nnod;
  *ndim = (*this)->ndim;
  *xyz  = (*this)->nodedata;
}

void
PyPF::Nodedata::setArray(int nnod, int ndim, double xyz[])
{
  int N = (*this)->nnod;
  int D = (*this)->ndim;
  double* _xyz = (*this)->nodedata;
  for (int n=0; n<N; n++)
    for (int d=0; d<D; d++)
      _xyz[n*D+d] = xyz[n*D+d];
}
