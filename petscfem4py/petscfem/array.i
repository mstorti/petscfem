// -*- c++ -*-
// $Id: array.i,v 1.1.2.6 2006/06/15 21:52:59 dalcinl Exp $

%include numpy.i

%inline %{
#define PyPF_INT    PyArray_INT
#define PyPF_FLOAT  PyArray_DOUBLE
%}

%define ARRAY_STDVEC_OUTPUT(VECARG, VALUE_T)
%typemap(in,numinputs=0) (VECARG) ($*ltype temp) "$1 = &temp;"
%typemap(argout) (VECARG)
{
  PyObject* o = ARRAY_NEW(&(*$1)[0], VALUE_T, 1, (*$1).size());
  ARRAY_arg_fail($symname, $argnum);
  %append_output(o);
}
%enddef
