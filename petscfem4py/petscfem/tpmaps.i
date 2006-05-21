// -*- c++ -*-
// $Id: tpmaps.i,v 1.1.2.2 2006/05/21 03:47:18 dalcinl Exp $


// pair of int values
%typemap(in, numinputs=0, noblock=1) (int*, int*)
($*1_ltype temp1=0, $*2_ltype temp2=0)
{ $1= &temp1; $2= &temp2; }
%typemap(argout, noblock=1) (int*, int*)
{ %append_output(Py_BuildValue("ii", *$1, *$2)) ;}


// PETSc objects
%define PETSC_OBJECT_TYPEMAP(Type)
typedef struct _p_##Type* Type;
%types(Type);
%typemap(typecheck) Type = SWIGTYPE*;
%typemap(arginit)   Type "$1 = PETSC_NULL;";
%typemap(check, noblock=1) Type
{ if ($1 == PETSC_NULL) %argument_nullref($type, $symname, $argnum); }
%typemap(in, noblock=1) Type (void  *argp = 0, int res = 0) {
  res = SWIG_ConvertPtr($input, &argp,$descriptor, $disown | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    PyObject* sobj = PyObject_GetAttr($input, SWIG_This());
    if (PyErr_Occurred()) PyErr_Clear();
    res = SWIG_ConvertPtr(sobj, &argp,$descriptor, $disown | %convertptr_flags);
    Py_XDECREF(sobj);
  }
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$*ltype", $symname, $argnum); 
  }
  $1 = %reinterpret_cast(argp, $ltype);
}
%typemap(freearg) Type "";

%enddef

