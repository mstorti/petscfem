// -*- c++ -*-
// $Id: array.i,v 1.1.2.4 2006/05/24 20:57:07 dalcinl Exp $

%include numpy.i

ARRAY_NUMTYPE(Int,   PyPF_INT,   PyArray_INT)
ARRAY_NUMTYPE(Float, PyPF_FLOAT, PyArray_DOUBLE)


%define ARRAY_STDVEC_OUTPUT(VECARG, VALUE_T)
%typemap(in,numinputs=0) (VECARG) ($*ltype temp) "$1 = &temp;"
%typemap(argout) (VECARG)
{
  PyObject* o = NULL;
  o = ARRAY_NEW(&(*$1)[0], VALUE_T, 1, (*$1).size());
  ARRAY_arg_fail($symname, $argnum);
  %append_output(o);
}
%enddef
