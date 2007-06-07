// -*- c++ -*-
// $Id$

%template() std::pair<int,int>;
%template() std::pair<char,int>;

%typemap(doc, name="array",type="int[]") 
(int rows, int cols, const int array[])
"$3_name: $3_type value";;
ARRAY_2D(int rows, int cols, const int array[], ARRAY_INPUT, PyPF_INT)
ARRAY_TYPECHECK_SEQUENCE((int rows, int cols, const int array[]), ARRAY_TYPECHECK_INT32)

%typemap(doc, name="array",type="double[]") 
(int rows, int cols, const double array[])
"$3_name: $3_type value";;
ARRAY_2D(int rows, int cols, const double array[], ARRAY_INPUT, PyPF_FLOAT)
ARRAY_TYPECHECK_SEQUENCE((int rows, int cols, const double array[]), ARRAY_TYPECHECK_FLOAT)

%include "DTable.h"

// PF4PY_NAMESPACE_BEGIN
// %feature("implicitconv") DTable<int>;
// PF4PY_NAMESPACE_END

// PF4PY_NAMESPACE_BEGIN
// %feature("implicitconv") DTable<double>;
// PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
%ignore DTable<int>::operator int*;
%ignore DTable<int>::operator const int*;
PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
%ignore DTable<double>::operator double*;
%ignore DTable<double>::operator const double*;
PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
%template(DTableI) DTable<int>;
%template(DTableS) DTable<double>;
PF4PY_NAMESPACE_END


// array interface <http://numeric.scipy.org/array_interface.html>
%header %{
namespace numpy{
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
%array_interface( PF4PY_NAMESPACE::DTable<int> );
%array_interface( PF4PY_NAMESPACE::DTable<double> );
