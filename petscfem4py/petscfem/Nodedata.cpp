// $Id: Nodedata.cpp,v 1.1.2.3 2006/03/06 16:56:04 rodrigop Exp $

#include "Nodedata.h"

#include <fem.h>

PyPF::Nodedata::~Nodedata()
{ }

PyPF::Nodedata::Nodedata()
{ }

PyPF::OptionTable*
PyPF::Nodedata::get_opt_table(bool create)
{
  OptionTable*& options = (*this)->options;
  if(options == NULL && create) options = new OptionTable;
  return options;
}

void
PyPF::Nodedata::getData(int* nnod, int* ndim, double* xyz[])
{
  /* test */
  if ((*this)->nodedata == NULL) throw Error("data not set");

  Nodedata::Base* nodedata = *this;
  *nnod = nodedata->nnod;
  *ndim = nodedata->ndim;
  *xyz  = nodedata->nodedata;
}

void
PyPF::Nodedata::setData(int nnod, int ndim, double xyz[])
{
  /* test */
  if (nnod < 1) throw Error("invalid number of nodes (nnod < 1)");
  if (ndim < 1) throw Error("invalid number of dimensions (ndim < 1)");
  if (ndim > 3) throw Error("invalid number of dimensions (ndim > 3)");
  
  Nodedata::Base*& nodedata = *this;
  if (nodedata == NULL) nodedata = new Nodedata::Base;

  if (nodedata->nodedata == NULL) {
    nodedata->nnod = nnod;
    nodedata->ndim = ndim;
    nodedata->nu   = ndim;
    nodedata->nodedata = new double[nnod*ndim];
  } else {
    if (nodedata->nnod != nnod) throw Error("invalid number of nodes");
    if (nodedata->ndim != ndim) throw Error("invalid number of dimensions");
  }
  memcpy(nodedata->nodedata, xyz, nnod*ndim*sizeof(double));
}

void
PyPF::Nodedata::getSize(int* nnod, int* ndim)
{
  /* test */
  if ((*this)->nodedata == NULL) throw Error("data not set");
  
  Nodedata::Base* nodedata = *this;
  *nnod = nodedata->nnod;
  *ndim = nodedata->ndim;
}

void
PyPF::Nodedata::view()
{
  Nodedata::Base* nodedata = *this;
  if (nodedata == NULL || nodedata->nodedata == NULL) return;

  int     nnod = nodedata->nnod;
  int     ndim = nodedata->ndim;
  double* data = nodedata->nodedata;
  for (int i=0; i<nnod; i++) {
    printf("%7d:", i);
    for (int j=0; j<ndim; j++) {
      printf(" %f", data[i*ndim+j]);
    }
    printf("\n");
  }
}
