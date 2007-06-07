// -*- c++ -*-
// $Id$

%runtime %{
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


// MPI communicator
// ================
%typemap(typecheck, match="in") MPI_Comm = SWIGTYPE;
%typemap(arginit, noblock=1)    MPI_Comm { $1 = MPI_COMM_NULL; }
%typemap(in, noblock=1)         MPI_Comm = SWIGTYPE;
%typemap(freearg, match="in")   MPI_Comm "";


// PETSc objects
// =============
%define %petsc_null_fmt(_type,_name,_argn) 
"null PETSc object " %argfail_fmt(_type, _name, _argn)
%enddef
%define %argument_petsc_null(type, name, argn) 
SWIG_exception_fail(SWIG_ValueError, %petsc_null_fmt(type, name, argn))
%enddef

%define PETSC_OBJECT_TYPEMAP(Type)
// ----
%types(Type);
typedef struct _p_##Type* Type;
// ----
%typemap(in) Type = SWIGTYPE;
%typemap(check, noblock=1) Type 
{ if ($1 == PETSC_NULL) %argument_petsc_null($type, $symname, $argnum); }
%typemap(typecheck, match="in") Type = SWIGTYPE;
%typemap(freearg, match="in") Type "";
// ----
%typemap(in, noblock=1) Type OPTIONAL (void *argp, int res = 0) {
  res = SWIG_ConvertPtr($input, &argp, $&descriptor, %convertptr_flags);
  if (!SWIG_IsOK(res)) %argument_fail(res, "$type", $symname, $argnum); 
  $1 = (!argp) ? PETSC_NULL : *(%reinterpret_cast(argp, $&ltype));
}
%typemap(check) Type OPTIONAL "";
// ----
%enddef
