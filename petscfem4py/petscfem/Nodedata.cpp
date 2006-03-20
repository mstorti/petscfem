// $Id: Nodedata.cpp,v 1.1.2.4 2006/03/20 16:06:00 rodrigop Exp $

#include "Nodedata.h"

#include <fem.h>


PYPF_NAMESPACE_BEGIN


OptionTable*
Nodedata::get_opt_table() const
{ 
  Nodedata::Base* nodedata = *this;
  if (nodedata->options == NULL) nodedata->options = new OptionTable;
  return nodedata->options;
}

Nodedata::~Nodedata()
{ 
  Nodedata::Base* nodedata = *this;
  PYPF_DELETE_VCTR(nodedata->nodedata);
  PYPF_DELETE_SCLR(nodedata->options);
  delete nodedata;
}

Nodedata::Nodedata()
  : Ptr(new Nodedata::Base), Object()
{ 
  Nodedata::Base* nodedata = *this;
  nodedata->nnod     = 0;
  nodedata->ndim     = 0;
  nodedata->nu       = 0;
  nodedata->nodedata = NULL;
  nodedata->options  = new OptionTable;
}

Nodedata::Nodedata(const Nodedata& _nodedata)
  : Ptr(new Nodedata::Base), Object(_nodedata)
{ 
  Nodedata::Base* nodedata = *this;

  nodedata->nnod     = 0;
  nodedata->ndim     = 0;
  nodedata->nu       = 0;
  nodedata->nodedata = NULL;
  nodedata->options  = new OptionTable;

  int nnod, ndim; double* xyz;
  _nodedata.getData(&nnod, &ndim, &xyz);
  this->setData(nnod, ndim, xyz);

  this->setOptions(_nodedata.getOptions());
}

Nodedata::Nodedata(Nodedata::Base* _nodedata)
  : Ptr(_nodedata), Object()
{  
  Nodedata::Base* nodedata = *this;
  if (nodedata->nodedata == NULL) {
    nodedata->nnod     = 0;
    nodedata->ndim     = 0;
    nodedata->nu       = 0;
  }
  if (nodedata->options == NULL) {
    nodedata->options = new OptionTable;
  }
}

Nodedata::Nodedata(int nnod, int ndim)
  : Ptr(new Nodedata::Base), Object()
{
  /* test */
  if (nnod < 1) throw Error("invalid number of nodes (nnod < 1)");
  if (ndim < 1) throw Error("invalid number of dimensions (ndim < 1)");
  if (ndim > 3) throw Error("invalid number of dimensions (ndim > 3)");
  Nodedata::Base* nodedata = *this;
  nodedata->nnod = nnod;
  nodedata->ndim = ndim;
  nodedata->nu   = ndim;
  nodedata->nodedata = new double[nnod*ndim];
  nodedata->options  = new OptionTable;
}  

Nodedata::Nodedata(int nnod, int ndim, double xyz[])
  : Ptr(new Nodedata::Base), Object()
{
  /* test */
  if (nnod < 1) throw Error("invalid number of nodes (nnod < 1)");
  if (ndim < 1) throw Error("invalid number of dimensions (ndim < 1)");
  if (ndim > 3) throw Error("invalid number of dimensions (ndim > 3)");
  /* */
  Nodedata::Base* nodedata = *this;
  nodedata->nnod = nnod;
  nodedata->ndim = ndim;
  nodedata->nu   = ndim;
  nodedata->nodedata = new double[nnod*ndim];
  memcpy(nodedata->nodedata, xyz, nnod*ndim*sizeof(double));
  nodedata->options  = new OptionTable;
}  


void
Nodedata::getData(int* nnod, int* ndim, double* xyz[]) const
{
  Nodedata::Base* nodedata = *this;
  *nnod = nodedata->nnod;
  *ndim = nodedata->ndim;
  *xyz  = nodedata->nodedata;
}

void
Nodedata::setData(int nnod, int ndim, double xyz[])
{
  /* test */
  if (nnod < 1) throw Error("invalid number of nodes (nnod < 1)");
  if (ndim < 1) throw Error("invalid number of dimensions (ndim < 1)");
  if (ndim > 3) throw Error("invalid number of dimensions (ndim > 3)");
  /* */
  Nodedata::Base* nodedata = *this;
  if (nodedata->nodedata == NULL) {
    nodedata->nnod = nnod;
    nodedata->ndim = ndim;
    nodedata->nu   = ndim;
    nodedata->nodedata = new double[nnod*ndim];
  } else {
    if (nodedata->nnod != nnod) throw Error("invalid number of nodes");
    if (nodedata->ndim != ndim) throw Error("invalid number of dimensions");
    if (nodedata->nu   != ndim) {
      nodedata->nu       = ndim;
      delete[] nodedata->nodedata;
      nodedata->nodedata = new double[nnod*ndim];
    }
  }
  memcpy(nodedata->nodedata, xyz, nnod*ndim*sizeof(double));
}

void
Nodedata::delData() 
{
  Nodedata::Base* nodedata = *this;
  PYPF_DELETE_VCTR(nodedata->nodedata);
  nodedata->nnod     = 0;
  nodedata->ndim     = 0;
  nodedata->nu       = 0;
  nodedata->nodedata = NULL;
}

void
Nodedata::getSize(int* nnod, int* ndim) const
{
  Nodedata::Base* nodedata = *this;
  *nnod = nodedata->nnod;
  *ndim = nodedata->ndim;
}

void
Nodedata::setUp()
{

}

void
Nodedata::clear() 
{
  Nodedata::Base* nodedata = *this;
  PYPF_DELETE_VCTR(nodedata->nodedata);
  PYPF_DELETE_SCLR(nodedata->options);
  nodedata->nnod     = 0;
  nodedata->ndim     = 0;
  nodedata->nu       = 0;
  nodedata->nodedata = NULL;
  nodedata->options  = new OptionTable;
}


void
Nodedata::view() const
{
  Nodedata::Base* nodedata = *this;

  int     nnod = nodedata->nnod;
  int     ndim = nodedata->ndim;
  double* data = nodedata->nodedata;

  for (int i=0; i<nnod; i++) {
    printf("Node %7d:", i);
    for (int j=0; j<ndim; j++) {
      printf(" %f", data[i*ndim+j]);
    }
    printf("\n");
  }

}


PYPF_NAMESPACE_END
