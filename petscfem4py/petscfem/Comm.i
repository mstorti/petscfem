// -*- c++ -*-
// $Id: Comm.i,v 1.1.2.1 2006/04/27 19:09:17 rodrigop Exp $


PYPF_NAMESPACE_BEGIN
%ignore Comm::operator=;
%ignore Comm::operator MPI_Comm*;
%ignore Comm::operator MPI_Comm&;
PYPF_NAMESPACE_END

%typemap(out, noblock=1) MPI_Comm&
{ %set_output(SWIG_NewPointerObj(%as_voidptr($1), $descriptor, %newpointer_flags)); }

PYPF_NAMESPACE_BEGIN
%extend Comm {
  MPI_Comm& getComm() { return *self; }
  void      setComm(MPI_Comm comm) { *self = comm; }
  %pythoncode {
  comm = property(getComm, setComm,  doc='MPI handle')
  size = property(getSize, doc='communicator size')
  rank = property(getRank, doc='communicator rank')
  }
}
PYPF_NAMESPACE_END

%include "Comm.h"
