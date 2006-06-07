// $Id: Nodeset.cpp,v 1.1.2.7 2006/06/07 16:27:21 dalcinl Exp $

#include "Nodeset.h"

#include <fem.h>

PYPF_NAMESPACE_BEGIN

void Nodeset::chk_sizes(int nnod, int ndim, int nval) const
{
  if (this->nodedata.empty()) {
    PYPF_ASSERT(nnod>=0, "invalid number of nodes (nnod<0)");
    PYPF_ASSERT(ndim>=1, "invalid number of dimensions, out of range (ndim<1)");
    PYPF_ASSERT(ndim<=3, "invalid number of dimensions, out of range (ndim>3)");
    PYPF_ASSERT(nval>=0, "invalid number of values, out of range (nval<0)");
  } else {
    PYPF_ASSERT(this->nnod==nnod, "cannot change number of nodes");
    PYPF_ASSERT(this->ndim==ndim, "cannot change number of dimensions");
    PYPF_ASSERT(this->nval==nval, "cannot change number of values");
  }
}

void Nodeset::set_sizes(int nnod, int ndim, int nval)
{ 
  this->nnod = nnod;
  this->ndim = ndim;
  this->nval = nval;
  this->nodedata.resize(nnod * (ndim + nval));
}

void Nodeset::set_array(const double array[])
{ 
  int     size = this->nodedata.size();
  double* data = &this->nodedata[0];
  if (size == 0) return;
  PYPF_ASSERT(array!=NULL, "null pointer to array data");
  memcpy(data, array, size*sizeof(double));
}

void Nodeset::touch() const {
  // get base pointer
  Nodeset::Base* nodedata = *this;
  // update options
  nodedata->options  = this->options;
  // update nodedata
  nodedata->nnod     = this->nnod;
  nodedata->ndim     = this->ndim;
  nodedata->nu       = this->ndim + this->nval;
  nodedata->nodedata = const_cast<double*>(&this->nodedata[0]);
}  

PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN

Nodeset::~Nodeset()
{ 
  Nodeset::Base* nodedata = *this;
  if (nodedata == NULL) return;
  nodedata->nodedata = NULL;
  nodedata->options  = NULL;
  delete nodedata;
}

Nodeset::Nodeset()
  : Handle(new Nodeset::Base),
    Object(),
    nnod(0), ndim(0), nval(0), nodedata()
{ 
  this->touch();
}

Nodeset::Nodeset(const Nodeset& nd)
  : Handle(new Nodeset::Base), 
    Object(nd),
    nnod(nd.nnod), ndim(nd.ndim), nval(nd.nval), nodedata(nd.nodedata)
{ 
  this->touch();
}

Nodeset::Nodeset(int nnod, int ndim, int nval) 
  : Handle(new Nodeset::Base), 
    Object(),
    nnod(0), ndim(0), nval(0), nodedata()
{
  this->chk_sizes(nnod, ndim, nval);
  this->set_sizes(nnod, ndim, nval);
  this->touch();
}

Nodeset::Nodeset(int nnod, int ndim, const double xnod[])
  : Handle(new Nodeset::Base),
    Object(),
    nnod(0), ndim(0), nval(0), nodedata()
{

  this->chk_sizes(nnod, ndim, 0);
  this->set_sizes(nnod, ndim, 0);
  this->set_array(xnod);
  this->touch();
}  

void
Nodeset::getSizes(int* nnod, int* ndim, int* nval) const
{
  if (nnod) *nnod = this->nnod;
  if (ndim) *ndim = this->ndim;
  if (nval) *nval = this->nval;
}

void 
Nodeset::setSizes(int nnod, int ndim, int nval) 
{
  this->chk_sizes(nnod, ndim, nval);
  this->set_sizes(nnod, ndim, nval);
  this->touch();
  
}

void
Nodeset::getArray(const double* array[]) const
{
  if (array) *array = &this->nodedata[0];
}

void
Nodeset::setArray(const double data[])
{

  this->set_array(data);
  this->touch();
}

int
Nodeset::getDim() const
{
  return this->ndim;
}

int
Nodeset::getSize() const
{
  return this->nnod;
}

void
Nodeset::getNode(int i, int* n, const double* node[]) const 
{
  int nnod, ndim, nval;  this->getSizes(&nnod, &ndim, &nval);
  const double* data;    this->getArray(&data);
  PYPF_ASSERT(i>=0 && i<nnod, "index out of range");
  if (n)    *n    = ndim + nval;
  if (data) *node = &data[i * (ndim + nval)];
}

void 
Nodeset::setNode(int i, int n, const double node[])
{
  int nnod, ndim, nval;  this->getSizes(&nnod, &ndim, &nval);
  const double* data;    this->getArray(&data);
  PYPF_ASSERT(i>=0 && i<nnod, "index out of range");
  PYPF_ASSERT(n==ndim || n==(ndim + nval) , "invalid number of dimensions/values");
  double* array = const_cast<double*>(&data[i * (ndim + nval)]);
  for (int k=0; k<n; k++) array[k] = node[k];
}

void
Nodeset::clear()
{
  this->options.clear();
  this->nodedata.clear();
  this->nnod = 0;
  this->touch();
}


void
Nodeset::view() const
{
  int nnod = this->nnod;
  int ndim = this->ndim;
  int nval = this->nval;
  int cols = ndim+nval;
  printf("Nodeset Object:\n");
  printf("  nnod=%d, ndim=%d, nval=%d\n", nnod, ndim, nval);
  int i=0;
  while (i<nnod) {
    printf("  ");
    printf(" %d -> ", i);
    for (int j=0; j<cols; j++) {
      if (j>0) {
	if (j == ndim) printf(" - "); 
	else           printf(", "); 
      }
      printf("%g", this->nodedata[i*cols+j]);
    }
    printf("\n");
    i++;
  }
}

// void 
// Nodeset::sync(int root) {
//   /* test */
//   int size;  MPI_Comm_size(this->comm, &size);
//   if (root < 0 || root >= size) 
//     throw Error("invalid root processor, out of range");
//   else if (size == 1) return; // nothing to do
//   // broadcast sizes
//   int ndsize[3];
//   ndsize[0] = this->nnod;
//   ndsize[1] = this->nval;
//   ndsize[2] = this->ndim;
//   MPI_Bcast(ndsize, 3, MPI_INT, root, this->comm);
//   int nnod = this->nnod = ndsize[0];
//   int nval = this->nval = ndsize[1];
//   int ndim = this->ndim = ndsize[2];
//   // broadcast data
//   int rank; MPI_Comm_rank(this->comm, &rank);
//   if (rank != root) this->nodedata.resize(nnod*nval);
//   double* buff = &this->nodedata[0]; int count = nnod*nval;
//   MPI_Bcast(buff, nnod*nval, MPI_DOUBLE, root, this->comm);

//   // base pointer
//   Nodeset::Base* nodedata = *this;
//   // options
//   nodedata->options  = this->options;
//   // nodedata
//   nodedata->nnod     = this->nnod;
//   nodedata->nu       = this->nval;
//   nodedata->ndim     = this->ndim;
//   nodedata->nodedata = &this->nodedata[0];
// }

PYPF_NAMESPACE_END
