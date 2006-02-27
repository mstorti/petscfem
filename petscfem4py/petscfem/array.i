// -*- c++ -*-

%include numpy.i

ARRAY_NUMTYPE(Int,   PyPF_INT,   PyArray_INT)
ARRAY_NUMTYPE(Float, PyPF_FLOAT, PyArray_DOUBLE)

%pythoncode %{
Int    = _petscfem.Int
Float  = _petscfem.Float
%}
