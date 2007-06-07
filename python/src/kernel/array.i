// -*- c++ -*-
// $Id$

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

%header %{
namespace numpy
{
  template<typename T> static char typechar() { return '\0'; }
  template<> static char typechar<int>()      { return  'i'; }
  template<> static char typechar<double>()   { return  'f'; }
}
%}

%define %array_interface(Class)
%extend Class { PyObject* __array_interface__; }
%{
#define %mangle(Class) ##_## __array_interface__ ## _get(_t) numpy::array_interface(_t)
#define %mangle(Class) ##_## __array_interface__ ## _set(_t, _val) \
        SWIG_exception_fail(SWIG_AttributeError, "read-only attribute")
%}
%enddef
