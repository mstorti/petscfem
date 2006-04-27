// $Id: Comm.cpp,v 1.1.2.1 2006/04/27 19:09:17 rodrigop Exp $


#include "Comm.h"
#include "Error.h"
#include "macros.h"

PYPF_NAMESPACE_BEGIN
 
Comm::Comm()
  : comm(MPI_COMM_NULL)
{ }

Comm::Comm(const Comm& comm)
  : comm(comm.comm)
{ }

Comm::Comm(MPI_Comm comm)
  : comm(comm)
{ }

int
Comm::getSize() const
{
  PYPF_ASSERT(this->comm!=MPI_COMM_NULL, "null communicator");
  int size;
  MPI_Comm_size(this->comm, &size);
  return size;
}

int
Comm::getRank() const
{
  PYPF_ASSERT(this->comm!=MPI_COMM_NULL, "null communicator");
  int rank;
  MPI_Comm_rank(this->comm, &rank);
  return rank;
}

bool 
Comm::operator==(const Comm& comm) const
{
  return this->comm == comm.comm;
}

bool 
Comm::operator!=(const Comm& comm) const
{
  return this->comm != comm.comm;
}

Comm& 
Comm::operator=(const Comm& comm)
{
  this->comm = comm.comm; return *this;
}

Comm& 
Comm::operator=(MPI_Comm comm) 
{
  this->comm = comm; return *this;
}


PYPF_NAMESPACE_END
