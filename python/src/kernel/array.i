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

// array interface <http://numeric.scipy.org/array_interface.html>
%header %{
namespace numpy {

  template<typename T> static char typechar() { return '\0'; }
  template<> static char typechar<int>()      { return  'i'; }
  template<> static char typechar<double>()   { return  'f'; }

  template<typename T>
  static PyObject* 
  array_interface(const PF4PY_NAMESPACE::DTable<T>* self)
  {
    const std::pair<int,int>& shape = self->getShape();
    const std::vector<T>&     data  = self->getArray();
    void* array = const_cast<T*>(&data[0]); 
    char endian = PyArray_NATIVE;
    char kind   = typechar<T>();
    int  elsize = sizeof(T);
    return Py_BuildValue("{sNsNsNsN}",
			 "shape",   Py_BuildValue("ii", shape.first, shape.second),
			 "typestr", PyString_FromFormat("%c%c%d", endian, kind, elsize),
			 "data",    Py_BuildValue("NO", PyLong_FromVoidPtr(array), Py_False),
			 "version", PyInt_FromLong(3));
  }
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
