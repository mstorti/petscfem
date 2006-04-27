// -*- c++ -*-

#ifndef PYPF_COMM_H
#define PYPF_COMM_H

#include <mpi.h>
#include "namespace.h"

PYPF_NAMESPACE_BEGIN

class Comm {
  
 protected:
  MPI_Comm comm;
  
 public:
  Comm();
  Comm(const Comm& comm);
  Comm(MPI_Comm comm);

  Comm& operator=(const Comm& comm);
  Comm& operator=(MPI_Comm comm);

  bool operator==(const Comm& comm) const;
  bool operator!=(const Comm& comm) const;
  
  int getSize() const;
  int getRank() const;
  
  inline operator MPI_Comm*() const 
  { return const_cast<MPI_Comm*>(&this->comm); }
  inline operator MPI_Comm&() const
  { return const_cast<MPI_Comm&>(this->comm); }

};

PYPF_NAMESPACE_END

#endif // PYPF_COMM_H
