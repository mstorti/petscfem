// -*- c++ -*-

PETSC_OBJECT_TYPEMAP(Vec)
PETSC_OBJECT_TYPEMAP(Mat)

//%typemap(check, noblock=1) Vec x0 "";
//%typemap(check, noblock=1) Vec x1 "";
%typemap(check, noblock=1) Mat J  "";

%include "NavierStokes.h"
