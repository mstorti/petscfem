// -*- c++ -*-
// $Id: tpmaps.i,v 1.1.2.6 2006/06/15 21:50:52 dalcinl Exp $

%header %{
SWIGINTERNINLINE PyObject* 
SWIG_getattr_this(PyObject* obj) {
  if (!obj) return NULL;
  obj = PyObject_GetAttr(obj, SWIG_This());
  if (!obj) PyErr_Clear();
  return obj;
}
SWIGINTERNINLINE int
SWIG_convert_ptr(PyObject *obj, void **ptr, swig_type_info *ty, int flags) {
  int res = SWIG_ConvertPtr(obj, ptr, ty, flags);
  if (!SWIG_IsOK(res)) {
    PyObject* _this = SWIG_getattr_this(obj);
    res = SWIG_ConvertPtr(_this, ptr, ty, flags);
    Py_XDECREF(_this);
  }
  return res;
}
#undef  SWIG_ConvertPtr
#define SWIG_ConvertPtr(obj, pptr, type, flags) \
        SWIG_convert_ptr(obj, pptr, type, flags)
%}

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


// MPI communicator
%typemap(typecheck, match="in") MPI_Comm = SWIGTYPE;
%typemap(arginit, noblock=1)    MPI_Comm { $1 = MPI_COMM_NULL; }
%typemap(in, noblock=1)         MPI_Comm = SWIGTYPE;
%typemap(freearg, match="in")   MPI_Comm "";


// PETSc objects
%define PETSC_OBJECT_TYPEMAP(Type)
typedef struct _p_##Type* Type;
%types(Type);
%typemap(typecheck, match="in") Type = SWIGTYPE*;
%typemap(arginit, noblock=1)    Type { $1 = PETSC_NULL; }
%typemap(in)                    Type = SWIGTYPE*;
%typemap(check, noblock=1)      Type
{ if ($1 == PETSC_NULL) %argument_nullref($type, $symname, $argnum); }
%typemap(freearg, match="in")   Type "";
%enddef
