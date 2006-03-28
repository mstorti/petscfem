// $Id: Nodedata.cpp,v 1.1.2.5 2006/03/28 22:13:25 rodrigop Exp $

#include "Nodedata.h"

#include <fem.h>


PYPF_NAMESPACE_BEGIN

Nodedata::~Nodedata()
{ 
  Nodedata::Base* nodedata = *this;
  nodedata->nodedata = NULL;
  nodedata->options  = NULL;
  delete nodedata;
}

Nodedata::Nodedata()
  : Handle(new Nodedata::Base), Object(),
    nnod(0), ndim(0), nodedata()
{ 
  // base pointer
  Nodedata::Base* nodedata = *this;
  // options
  nodedata->options  = this->options;
  // nodedata
  nodedata->nnod     = this->nnod;
  nodedata->ndim     = this->ndim;
  nodedata->nu       = this->ndim;
  nodedata->nodedata = &this->nodedata[0];
}

Nodedata::Nodedata(const Nodedata& nd)
  : Handle(new Nodedata::Base), Object(nd),
    nnod(nd.nnod), ndim(nd.ndim), nodedata(nd.nodedata)
{ 
  Nodedata::Base* nodedata = *this;
  // options
  nodedata->options  = this->options;
  // nodedata
  nodedata->nnod     = this->nnod;
  nodedata->ndim     = this->ndim;
  nodedata->nu       = this->ndim;
  nodedata->nodedata = &this->nodedata[0];
}

Nodedata::Nodedata(Nodedata::Base* base)
  : Handle(base), Object(),
    nnod(base->nnod), ndim(base->ndim), nodedata(nnod*ndim)
{  
  // nodedata
  if (base->nodedata != NULL) {
    memcpy(&this->nodedata[0], base->nodedata,
	   this->nnod*this->ndim*sizeof(double));
    delete[] base->nodedata;
  }
  // base pointer
  Nodedata::Base* nodedata = *this;
  // nodedata
  nodedata->nnod     = this->nnod;
  nodedata->ndim     = this->ndim;
  nodedata->nu       = this->ndim;
  nodedata->nodedata = &this->nodedata[0];
  // options
  if (nodedata->options == NULL) 
    nodedata->options = this->options;
  else 
    this->options = nodedata->options;
}

Nodedata::Nodedata(int nnod, int ndim, const double xyz[])
  : Handle(new Nodedata::Base), Object(),
    nnod(nnod), ndim(ndim), nodedata(nnod*ndim)
{
  if (nnod*ndim !=0) {
    if (nnod < 1) throw Error("invalid number of nodes (nnod < 1)");
    if (ndim < 1) throw Error("invalid number of dimensions (ndim < 1)");
    if (ndim > 3) throw Error("invalid number of dimensions (ndim > 3)");
  } else {
    this->nnod = this->ndim = nnod = ndim = 0;
  }
  memcpy(&this->nodedata[0], xyz, nnod*ndim*sizeof(double));
  
  // base pointer
  Nodedata::Base* nodedata = *this;
  // options
  nodedata->options  = this->options;
  // nodedata
  nodedata->nnod     = this->nnod;
  nodedata->ndim     = this->ndim;
  nodedata->nu       = this->ndim;
  nodedata->nodedata = &this->nodedata[0];
}  


void
Nodedata::getSize(int* nnod, int* ndim) const
{
  if (nnod) *nnod = this->nnod;
  if (ndim) *ndim = this->ndim;
}

void
Nodedata::getData(int* nnod, int* ndim, double* xyz[]) const
{
  if (nnod) *nnod = this->nnod;
  if (ndim) *ndim = this->ndim;
  if (xyz)  *xyz  = const_cast<double*>(&this->nodedata[0]);
}

void
Nodedata::setData(int nnod, int ndim, const double xyz[])
{
  /* test */
  if (nnod*ndim !=0) {
    if (nnod < 1) throw Error("invalid number of nodes (nnod < 1)");
    if (ndim < 1) throw Error("invalid number of dimensions (ndim < 1)");
    if (ndim > 3) throw Error("invalid number of dimensions (ndim > 3)");
  } else  nnod = ndim = 0;
  
  if (this->nnod*this->ndim == 0) {
    this->nnod  = nnod;
    this->ndim  = ndim;
  } else {
    if (this->nnod != nnod) throw Error("invalid number of nodes");
    if (this->ndim != ndim) throw Error("invalid number of dimensions");
  }
  this->nodedata.resize(nnod*ndim);
  memcpy(&this->nodedata[0], xyz, nnod*ndim*sizeof(double));
  
  // base pointer
  Nodedata::Base* nodedata = *this;
  // nodedata
  nodedata->nnod     = this->nnod;
  nodedata->ndim     = this->ndim;
  nodedata->nu       = this->ndim;
  nodedata->nodedata = &this->nodedata[0];
}


std::vector<double> 
Nodedata::getNode(int n) const 
{
  if (n<0 || n>=this->nnod) throw Error("index out of range");
  const double* beg = &this->nodedata[n*this->ndim];
  const double* end = beg + this->ndim;
  return std::vector<double>(beg, end);
}

void 
Nodedata::setNode(int n, const std::vector<double>& node)
{
  if (n<0 || n>=this->nnod) throw Error("index out of range");
  if (node.size() != this->ndim) throw Error("invalid number of dimensions");
  double* data = &this->nodedata[n*this->ndim];
  for (int i=0; i<this->ndim; i++) data[i] = node[i];
}

void
Nodedata::setUp()
{ }

void
Nodedata::clear() 
{
  this->options.clear();
  this->nodedata.clear();
  this->nnod = 0;
  this->ndim = 0;
  
  // base pointer
  Nodedata::Base* nodedata = *this;
  // options
  nodedata->options  = this->options;
  // nodedata
  nodedata->nnod     = this->nnod;
  nodedata->ndim     = this->ndim;
  nodedata->nu       = this->ndim;
  nodedata->nodedata = &this->nodedata[0];
}


void
Nodedata::view() const
{
  int nnod = this->nnod;
  int ndim = this->ndim;
  for (int i=0; i<nnod; i++) {
    printf("Node %7d:", i);
    for (int j=0; j<ndim; j++) {
      printf(" %f", this->nodedata[i*ndim+j]);
    }
    printf("\n");
  }
}

void 
Nodedata::sync(int root) {
  /* test */
  int size;  MPI_Comm_size(this->comm, &size);
  if (root < 0 || root >= size) 
    throw Error("invalid root processor, out of range");
  else if (size == 1) return; // nothing to do
  // broadcast sizes
  int ndsize[2];
  ndsize[0] = this->nnod;
  ndsize[1] = this->ndim;
  MPI_Bcast(ndsize, 2, MPI_INT, root, this->comm);
  int nnod = this->nnod = ndsize[0];
  int ndim = this->ndim = ndsize[1];
  // broadcast data
  int rank; MPI_Comm_rank(this->comm, &rank);
  if (rank != root) this->nodedata.resize(nnod*ndim);
  double* buff = &this->nodedata[0]; int count = nnod*ndim;
  MPI_Bcast(buff, count, MPI_DOUBLE, root, this->comm);

  // base pointer
  Nodedata::Base* nodedata = *this;
  // options
  nodedata->options  = this->options;
  // nodedata
  nodedata->nnod     = this->nnod;
  nodedata->ndim     = this->ndim;
  nodedata->nu       = this->ndim;
  nodedata->nodedata = &this->nodedata[0];
}


PYPF_NAMESPACE_END
