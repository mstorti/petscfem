// $Id$

#ifndef PF4PY_COMM_H
#define PF4PY_COMM_H

#include <mpi.h>
#include "namespace.h"
#include "Error.h"

PF4PY_NAMESPACE_BEGIN

class Comm
#if defined(SWIG)
  : public MPI_Comm  // vile hack
#endif
{
private:
  MPI_Comm comm;

protected:
  inline bool eq(MPI_Comm other_comm) const {
    if (this->comm == MPI_COMM_NULL && other_comm == MPI_COMM_NULL) return true;
    if (this->comm == MPI_COMM_NULL || other_comm == MPI_COMM_NULL) return false;
    int result; MPI_Comm_compare(this->comm, other_comm, &result);
    if (result != MPI_IDENT && result != MPI_CONGRUENT) return false;
    return true;
  }
  inline void check() const { 
    if (this->comm == MPI_COMM_NULL)
      throw Error("Comm: null communicator"); 
  }
  
public:
  inline ~Comm() { }
  inline Comm()                 : comm(MPI_COMM_NULL) { }
  inline Comm(const Comm& comm) : comm(comm.comm)     { }
  inline Comm(MPI_Comm comm)    : comm(comm)          { }

  inline Comm& operator=(MPI_Comm comm) { this->comm = comm; return *this; }
  inline operator MPI_Comm() const { return this->comm; }

  inline bool operator==(MPI_Comm comm) const { return this->eq(comm); }
  inline bool operator!=(MPI_Comm comm) const { return !this->eq(comm); }
  
  inline int getSize() const { this->check(); int s; MPI_Comm_size(*this, &s); return s; }
  inline int getRank() const { this->check(); int r; MPI_Comm_rank(*this, &r); return r; }
};

PF4PY_NAMESPACE_END

#endif // PF4PY_COMM_H

// Local Variables:
// mode: C++
// End:
