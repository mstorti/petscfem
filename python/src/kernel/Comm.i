// -*- c++ -*-
// $Id$

//PF4PY_NAMESPACE_BEGIN
//%feature("implicitconv") Comm;
//PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
%ignore Comm::operator MPI_Comm;
%ignore Comm::operator MPI_Comm*;
PF4PY_NAMESPACE_END

// %typemap(out, noblock=1) MPI_Comm&
// { %set_output(SWIG_NewPointerObj(%as_voidptr($1), $descriptor, %newpointer_flags)); }

// PF4PY_NAMESPACE_BEGIN
// %extend Comm {
//   MPI_Comm& getComm() { return *$self; }
//   %pythoncode {
//   comm = property(getComm, doc='MPI handle (MPI_Comm)')
//   size = property(getSize, doc='communicator size')
//   rank = property(getRank, doc='communicator rank')
//   }
// }
// PF4PY_NAMESPACE_END

%include property.i
%property(PETScFEM::Comm, int, size, "communicator size", self->getSize() );
%property(PETScFEM::Comm, int, rank, "processor rank",    self->getRank() );

%ignore MPI_Comm;
struct MPI_Comm { };

%include "Comm.h"
