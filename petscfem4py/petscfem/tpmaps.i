// -*- c++ -*-
// $Id: tpmaps.i,v 1.1.2.5 2006/06/07 16:25:32 dalcinl Exp $

// pair of int values
%typemap(in, numinputs=0, noblock=1) (int*, int*)
($*1_ltype temp1=0, $*2_ltype temp2=0)
{ $1= &temp1; $2= &temp2; }
%typemap(argout, noblock=1) (int*, int*)
{ %append_output(Py_BuildValue("ii", *$1, *$2)) ;}

// triad of int values
%typemap(in, numinputs=0, noblock=1) (int*, int*, int*)
($*1_ltype temp1=0, $*2_ltype temp2=0, $*3_ltype temp3=0)
{ $1= &temp1; $2= &temp2; $3= &temp3; }
%typemap(argout, noblock=1) (int*, int*, int*)
{ %append_output(Py_BuildValue("iii", *$1, *$2, *$3)) ;}


%define %get_swig_this(obj, argp)
do { 
  *(argp) = (obj)? PyObject_GetAttr((obj), SWIG_This()) : 0;
  if (PyErr_Occurred()) PyErr_Clear();
} while(0)
%enddef


%typemap(in) MPI_Type (void *argp, int res = 0) {
  res = SWIG_ConvertPtr($input, &argp, $&descriptor, %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    PyObject* _swig_this; %get_swig_this($input, &_swig_this);
    res = SWIG_ConvertPtr(_swig_this, &argp, $&descriptor, %convertptr_flags);
    Py_XDECREF(_swig_this);
    if (!SWIG_IsOK(res)) { %argument_fail(res, "$type", $symname, $argnum); }  
  }
  if (argp) { $1 = *(%reinterpret_cast(argp, $&ltype));} 
  /*else    {  %argument_nullref("$type", $symname, $argnum); } */
}


// MPI communicator
%typemap(arginit, noblock=1) MPI_Comm { $1 = MPI_COMM_NULL; }
%typemap(typecheck) MPI_Comm = SWIGTYPE;
%typemap(in) MPI_Comm = MPI_Type;
%typemap(freearg) MPI_Comm "";



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
    PyObject* _swig_this; %get_swig_this($input, &_swig_this);
    res = SWIG_ConvertPtr(_swig_this, &argp,$descriptor, $disown | %convertptr_flags);
    Py_XDECREF(_swig_this);
    if (!SWIG_IsOK(res)) { %argument_fail(res, "$*ltype", $symname, $argnum);}
  }
  $1 = %reinterpret_cast(argp, $ltype);
}
%typemap(freearg) Type "";

%enddef

