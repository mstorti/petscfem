// -*- c++ -*-
// $Id$

%template() std::pair<std::string,int>;
%template() std::vector<std::pair<std::string,int> >;

%include "PTable.h"

PF4PY_NAMESPACE_BEGIN
%ignore PTable<int>::operator int*;
%ignore PTable<int>::operator const int*;
PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
%ignore PTable<double>::operator double*;
%ignore PTable<double>::operator const double*;
PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
%template(PTableI) PTable<int>;
%template(PTableS) PTable<double>;
PF4PY_NAMESPACE_END


// array interface <http://numeric.scipy.org/array_interface.html>
%header %{
namespace numpy{
  template<typename T>
  static PyObject* 
  array_interface(const PF4PY_NAMESPACE::PTable<T>* self)
  {
    typedef std::pair<std::string,int> FieldEntry;
    typedef std::vector<FieldEntry>    FieldList;
    //const FieldList&          fields = self->getFields();
    const std::pair<int,int>& shape  = self->getShape();
    const std::vector<T>&     data   = self->getArray();
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
%array_interface( PF4PY_NAMESPACE::PTable<int> );
%array_interface( PF4PY_NAMESPACE::PTable<double> );


