// -*- c++ -*-

#ifndef PYPF_COMM_H
#define PYPF_COMM_H

#include <mpi.h>
#include "namespace.h"
#include "Error.h"

#if defined(SWIG)
%ignore MPI_Comm;
class MPI_Comm { };
#endif

PYPF_NAMESPACE_BEGIN

class Comm
#if defined(SWIG)
  : public MPI_Comm
#endif
{
private:
  inline bool comp(MPI_Comm other_comm) const {
    if (this->comm == MPI_COMM_NULL && other_comm == MPI_COMM_NULL) return true;
    if (this->comm == MPI_COMM_NULL || other_comm == MPI_COMM_NULL) return false;
    int result; MPI_Comm_compare(this->comm, other_comm, &result);
    if (result == MPI_IDENT || result == MPI_CONGRUENT) return true;
    else return false;
  }
  inline void check() const 
  { if (this->comm == MPI_COMM_NULL) throw Error("null communicator"); }

protected:
  MPI_Comm comm;
  
public:
  ~Comm() { this->comm = MPI_COMM_NULL; };
  Comm()                 : comm(MPI_COMM_NULL) { }
  Comm(const Comm& comm) : comm(comm.comm)     { }
  Comm(MPI_Comm comm)    : comm(comm)          { }

  Comm& operator=(MPI_Comm comm)    { this->comm = comm; return *this; }
  inline operator MPI_Comm*() const { return const_cast<MPI_Comm*>(&this->comm); }
  inline operator MPI_Comm&() const { return const_cast<MPI_Comm&>(this->comm); }

  bool operator==(MPI_Comm comm) const { return this->comp(comm); };
  bool operator!=(MPI_Comm comm) const { return this->comp(comm); };
  
  int getSize() const { this->check(); int s; MPI_Comm_size(*this, &s); return s; }
  int getRank() const { this->check(); int r; MPI_Comm_rank(*this, &r); return r; }
};

PYPF_NAMESPACE_END

#endif // PYPF_COMM_H
