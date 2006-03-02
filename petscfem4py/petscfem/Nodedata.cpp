// $Id: Nodedata.cpp,v 1.1.2.2 2006/03/02 21:37:12 rodrigop Exp $

#include "Nodedata.h"

#include <fem.h>

PyPF::Nodedata::Nodedata() : Ptr(new ::Nodedata)
{ 
#if 0
  (*this)->nnod = 0;
  (*this)->ndim = 0;
  (*this)->nu   = 0;
  (*this)->nodedata = NULL;
  (*this)->options  = NULL;
#else
  
#endif
}

PyPF::Nodedata::~Nodedata()
{ 
#if 0
  if ((*this)->options != NULL) {
    delete (*this)->options; (*this)->options = NULL;
  }
  if (this->ptr != NULL) {
    delete this->ptr;
  }
#else

#endif
}


void
PyPF::Nodedata::getData(int* nnod, int* ndim, double* xyz[])
{
  if ((*this)->nodedata == NULL) {
    throw Error("node data not set");
  }
  *nnod = (*this)->nnod;
  *ndim = (*this)->ndim;
  *xyz  = (*this)->nodedata;
}

void
PyPF::Nodedata::setData(int nnod, int ndim, double xyz[])
{
  if ((*this)->nodedata == NULL) {
    throw Error("null node data");
    if (ndim < 1 || ndim > 3) throw Error("invalid number of dimensions");
    (*this)->nnod = nnod;
    (*this)->ndim = ndim;
    (*this)->nu   = ndim;
    (*this)->nodedata = new double[nnod*ndim];
    this->setData(nnod, ndim, xyz);
  } else {
    int N = (*this)->nnod; 
    int D = (*this)->ndim;
    if (N != nnod) throw Error("invalid number of nodes");
    if (D != ndim) throw Error("invalid number of dimensions");
    double* _xyz = (*this)->nodedata;
    for (int n=0; n<N; n++)
      for (int d=0; d<D; d++)
	_xyz[n*D+d] = xyz[n*D+d];
  }
}

void
PyPF::Nodedata::getSize(int* nnod, int* ndim) 
{
  if ((*this)->nodedata == NULL) {
    throw Error("node data not set");
  }
  *nnod = (*this)->nnod;
  *ndim = (*this)->ndim;
}


std::string
PyPF::Nodedata::getOption(const std::string& key)
{
  if ((*this)->options == NULL) {
    throw Error("null option table");
  }
  const char* value = NULL;
  (*this)->options->get_entry(key.c_str(), value);
  if (value == NULL) throw Error("option not found");
  return value;
}

void 
PyPF::Nodedata::setOption(const std::string& key,
			  const std::string& value)
{
  if ((*this)->options == NULL) {
    throw Error("null option table");
    (*this)->options = new TextHashTable;
  }
  (*this)->options->add_entry(key.c_str(), value.c_str());
}

void
PyPF::Nodedata::view()
{
  if ((*this) == NULL) return;
  if ((*this)->nodedata == NULL) return;
  int     nnod = (*this)->nnod;
  int     ndim = (*this)->ndim;
  double* data = (*this)->nodedata;
  for (int i=0; i<nnod; i++) {
    printf("%7d:", i);
    for (int j=0; j<ndim; j++) {
      printf(" %f", data[i*ndim+j]);
    }
    printf("\n");
  }
}
