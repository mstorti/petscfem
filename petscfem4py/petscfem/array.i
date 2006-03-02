// -*- c++ -*-
// $Id: array.i,v 1.1.2.2 2006/03/02 21:37:12 rodrigop Exp $

%include numpy.i

ARRAY_NUMTYPE(Int,   PyPF_INT,   PyArray_INT)
ARRAY_NUMTYPE(Float, PyPF_FLOAT, PyArray_DOUBLE)

%pythoncode %{
## Numeric Types
## -------------
Int    = _petscfem.Int
Float  = _petscfem.Float
%}
