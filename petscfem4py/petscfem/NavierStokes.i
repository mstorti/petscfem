// -*- c++ -*-

PETSC_OBJECT_TYPEMAP(Vec)
PETSC_OBJECT_TYPEMAP(Mat)

PYPF_NAMESPACE_BEGIN

%newobject NavierStokes::read;

PYPF_NAMESPACE_END

%include "NavierStokes.h"
