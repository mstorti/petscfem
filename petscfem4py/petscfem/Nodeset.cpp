// $Id: Nodeset.cpp,v 1.1.2.3 2006/05/24 21:10:46 dalcinl Exp $

#include "Nodeset.h"

#include <fem.h>


PYPF_NAMESPACE_BEGIN

Nodeset::~Nodeset()
{ 
  Nodeset::Base* nodedata = *this;
  nodedata->nodedata = NULL;
  nodedata->options  = NULL;
  delete nodedata;
}

Nodeset::Nodeset()
  : Handle(new Nodeset::Base), 
    Object(),
    ndim(0), nnod(0), nval(0), nodedata()
{ 
  // base pointer
  Nodeset::Base* nodedata = *this;
  // options
  nodedata->options  = this->options;
  // nodedata
  nodedata->ndim     = this->ndim;
  nodedata->nnod     = this->nnod;
  nodedata->nu       = this->nval;
  nodedata->nodedata = &this->nodedata[0];
}

Nodeset::Nodeset(const Nodeset& nd)
  : Handle(new Nodeset::Base), 
    Object(nd),
    ndim(nd.ndim), nnod(nd.nnod), nval(nd.nval), nodedata(nd.nodedata)
{ 
  // base pointer
  Nodeset::Base* nodedata = *this;
  // options
  nodedata->options  = this->options;
  // nodedata
  nodedata->nnod     = this->nnod;
  nodedata->nu       = this->nval;
  nodedata->ndim     = this->ndim;
  nodedata->nodedata = &this->nodedata[0];
}

Nodeset::Nodeset(int ndim, int nnod, int nval, const double data[])
  : Handle(new Nodeset::Base),
    Object(),
    ndim(0), nnod(0), nval(0), nodedata()
{
  
  this->setDim(ndim);
  this->setData(nnod, nval, data);

  // base pointer
  Nodeset::Base* nodedata = *this;
  // options
  nodedata->options  = this->options;
  // nodedata
  nodedata->nnod     = this->nnod;
  nodedata->nu       = this->nval;
  nodedata->ndim     = this->ndim;
  nodedata->nodedata = &this->nodedata[0];
}  

void
Nodeset::setDim(int ndim)
{
  PYPF_ASSERT(this->ndim==0 || this->ndim==ndim,
	      "cannot change number of dimensions");
  PYPF_ASSERT(ndim>=1,
	      "invalid number of dimensions, out of range (ndim<1)");
  PYPF_ASSERT(ndim<=3, 
	      "invalid number of dimensions, out of range (ndim>3)");
  PYPF_ASSERT(this->nval==0 || ndim<=this->nval,
	      "invalid number of dimensions for data size (nval<ndim)");
  this->ndim = ndim;
}

int
Nodeset::getDim() const
{
  return this->ndim;
}

void
Nodeset::getData(int* nnod, int* nval, const double* data[]) const
{
  if (nnod) *nnod = this->nnod;
  if (nval) *nval = this->nval;
  if (data) *data = &this->nodedata[0];
}

void
Nodeset::setData(int nnod, int nval, const double data[])
{
  // test data
  if (nnod*nval != 0) {
    int ndim = this->ndim;
    PYPF_ASSERT(nnod>=1,    "invalid number of nodes (nnod<1)");
    PYPF_ASSERT(nval>=1,    "invalid number of values (nval<1)");
    PYPF_ASSERT(data!=NULL, "null pointer to data array");
    PYPF_ASSERT(nval>=ndim, "invalid number of values (nval<ndim)");
  } else nnod = nval = 0;
  if (this->nnod*this->nval != 0) {
    PYPF_ASSERT(this->nnod==nnod, "cannot change data size (nnod)");
    PYPF_ASSERT(this->nval==nval, "cannot change data size (nval)");
  }
  // copy data
  this->nnod  = nnod;
  this->nval  = nval;
  this->nodedata.resize(nnod*nval);
  memcpy(&this->nodedata[0], data, nnod*nval*sizeof(double));
  
  // base pointer
  Nodeset::Base* nodedata = *this;
  // options
  nodedata->options  = this->options;
  // nodedata
  nodedata->ndim     = this->ndim;
  nodedata->nnod     = this->nnod;
  nodedata->nu       = this->nval;
  nodedata->nodedata = &this->nodedata[0];
}

void
Nodeset::getDataSize(int* nnod, int* nval) const
{
  if (nnod) *nnod = this->nnod;
  if (nval) *nval = this->nval;
}

void
Nodeset::getNode(int i, int* n, const double* node[]) const 
{
  int nnod, nval; const double* data;
  this->getData(&nnod, &nval, &data);
  PYPF_ASSERT(i>=0 && i<nnod, "index out of range");
  if (n)    *n    = nval;
  if (data) *node = &data[i*nval];
}

void 
Nodeset::setNode(int i, int n, const double node[])
{
  int nnod, nval; const double* nodedata;
  this->getData(&nnod, &nval, &nodedata);
  PYPF_ASSERT(i>=0 && i<nnod, "index out of range");
  PYPF_ASSERT(n==nval || n==this->ndim , "invalid number of values");
  double* data = const_cast<double*>(&nodedata[i*nval]);
  for (int k=0; k<n; k++) data[k] = node[k];
}

int
Nodeset::getSize() const
{
  int nnod;
  this->getDataSize(&nnod, NULL);
  return nnod;
}

void
Nodeset::clear() 
{
  this->options.clear();
  this->nodedata.clear();
  this->nnod = 0;
  this->nval = 0;
  this->ndim = 0;
  
  // base pointer
  Nodeset::Base* nodedata = *this;
  // options
  nodedata->options  = this->options;
  // nodedata
  nodedata->ndim     = this->ndim;
  nodedata->nnod     = this->nnod;
  nodedata->nu       = this->nval;
  nodedata->nodedata = &this->nodedata[0];
}


void
Nodeset::view() const
{
  int nnod = this->nnod;
  int nval = this->nval;
  int ndim = this->ndim;
  printf("Nodeset Object:\n");
  printf("  size: nnod=%d, ndim=%d\n", nnod, ndim);
  printf("  data:");
  int i=0;
  while (i<nnod) {
    if (i>0) printf("       "); 
    printf(" %d -> ", i);
    for (int j=0; j<nval; j++) {
      if (j>0) printf(", "); 
      printf("%g", this->nodedata[i*nval+j]);
    }
    printf("\n");
    i++;
  }
}

void 
Nodeset::sync(int root) {
  /* test */
  int size;  MPI_Comm_size(this->comm, &size);
  if (root < 0 || root >= size) 
    throw Error("invalid root processor, out of range");
  else if (size == 1) return; // nothing to do
  // broadcast sizes
  int ndsize[3];
  ndsize[0] = this->nnod;
  ndsize[1] = this->nval;
  ndsize[2] = this->ndim;
  MPI_Bcast(ndsize, 3, MPI_INT, root, this->comm);
  int nnod = this->nnod = ndsize[0];
  int nval = this->nval = ndsize[1];
  int ndim = this->ndim = ndsize[2];
  // broadcast data
  int rank; MPI_Comm_rank(this->comm, &rank);
  if (rank != root) this->nodedata.resize(nnod*nval);
  double* buff = &this->nodedata[0]; int count = nnod*nval;
  MPI_Bcast(buff, nnod*nval, MPI_DOUBLE, root, this->comm);

  // base pointer
  Nodeset::Base* nodedata = *this;
  // options
  nodedata->options  = this->options;
  // nodedata
  nodedata->nnod     = this->nnod;
  nodedata->nu       = this->nval;
  nodedata->ndim     = this->ndim;
  nodedata->nodedata = &this->nodedata[0];
}

PYPF_NAMESPACE_END
