// $Id: Comm.cpp,v 1.1.2.2 2006/06/05 20:44:58 dalcinl Exp $


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

static bool comm_compare(MPI_Comm comm1, MPI_Comm comm2)
{
  if (comm1 == MPI_COMM_NULL && comm2 == MPI_COMM_NULL) return true;
  if (comm1 == MPI_COMM_NULL || comm2 == MPI_COMM_NULL) return false;
  int result;
  MPI_Comm_compare(comm1, comm2, &result);
  if (result == MPI_IDENT || result == MPI_CONGRUENT) return true;
  else return false;
}

bool 
Comm::operator==(const Comm& comm) const
{
  return comm_compare(this->comm, comm.comm);
}

bool 
Comm::operator==(MPI_Comm comm) const
{
  return comm_compare(this->comm, comm);
}

bool 
Comm::operator!=(const Comm& comm) const
{
  return !comm_compare(this->comm, comm.comm);
}

bool 
Comm::operator!=(MPI_Comm comm) const
{
  return !comm_compare(this->comm, comm);
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
