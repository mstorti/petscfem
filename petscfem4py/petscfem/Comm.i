// -*- c++ -*-
// $Id: Comm.i,v 1.1.2.3 2006/06/09 14:23:18 dalcinl Exp $

%header %{
#include "Comm.h"
%}

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
  %pythoncode {
  comm = property(getComm, doc='MPI handle (MPI_Comm)')
  size = property(getSize, doc='communicator size')
  rank = property(getRank, doc='communicator rank')
  }
}
PYPF_NAMESPACE_END

%include "Comm.h"
