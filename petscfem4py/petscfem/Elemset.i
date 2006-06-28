// -*- c++ -*-
// $Id: Elemset.i,v 1.1.2.8 2006/06/28 20:14:10 dalcinl Exp $


%include Object.i

%typemap(doc,name="elem",type="int[]")
  (int n, const int elem[])
  "elem: int[] value";
ARRAY_FLAT(int n, const int elem[],
	   ARRAY_INPUT, PyPF_INT)
ARRAY_TYPECHECK_SEQUENCE((int n, const int elem[]), 
			 ARRAY_TYPECHECK_INT32)
ARRAY_1D_NEW(int* n, int* elem[], PyPF_INT)


%typemap(doc,name="elemdata",type="int[]")
  (int nelem, int nel, const int icone[])
  "elemdata: int[] value";
ARRAY_2D(int nelem, int nel, const int icone[],
	 ARRAY_INPUT, PyPF_INT)
ARRAY_TYPECHECK_SEQUENCE((int nelem, int nel, const int icone[]),
			 ARRAY_TYPECHECK_INT32)
ARRAY_2D_NEW(int* nelem, int* nel, const int* icone[], PyPF_INT)


ARRAY_1D_NEW(int* n, int* part[], PyPF_INT)

%typemap(doc, type="list<Elemset>")
  std::vector<PYPF_NAMESPACE::Elemset*>&,
  const std::vector<PYPF_NAMESPACE::Elemset*>&
  "$1_name: list of Elemset instances";

%template() std::vector<PYPF_NAMESPACE::Elemset*>;
%typemap(check, noblock=1) 
  std::vector<PYPF_NAMESPACE::Elemset*>&,
  const std::vector<PYPF_NAMESPACE::Elemset*>& {
  for ($*1_ltype::const_iterator p = $1->begin(); p != $1->end(); p++) {
    if (*p == NULL) %argument_nullref("Elemset", $symname, $argnum);
  }
}


%apply (int*, int*) { (int* nelem, int* nel) };


PYPF_NAMESPACE_BEGIN
%pythonappend Elemset::getData %{val-=1%}
%pythonappend Elemset::getElem %{val-=1%}
PYPF_NAMESPACE_END

PYPF_NAMESPACE_BEGIN
%pythonappend Elemset::getPart %{val-=1%}
PYPF_NAMESPACE_END

PYPF_NAMESPACE_BEGIN
%extend Elemset {
  int __len__() 
    { return self->getSize(); }
  void __getitem__(int i, int* n, const int* elem[]) 
    { self->getElem(i, n, elem); }
  void __setitem__(int i, int  n, const int  elem[]) 
    { self->setElem(i, n, elem); }
  %pythoncode {
  def __iter__(self):
      """__iter__(self) -> iterator"""
      for i in xrange(len(self)):
          yield self[i]
  }
}
PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN
%extend Elemset {
  %pythoncode {
  type = property(getType,          doc='elemset type')
  name = property(getName, setName, doc='elemset name')
  size = property(getSize,          doc='number of elements')
  data = property(getData, setData, doc='elemset data')
  }
}
PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN
%extend Elemset {
  PyObject* __array_interface__ () {
    int nelem, nel;
    const int* icone; 
    self->getData(&nelem, &nel, &icone);
    char endian = PyArray_NATIVE;
    char kind   = PyArray_INTLTR;
    int  elsize = sizeof(int);
    return Py_BuildValue("{sNsNsNsN}",
			 "shape",   Py_BuildValue("ii", nelem, nel),
			 "typestr", PyString_FromFormat("%c%c%d", endian, kind, elsize),
			 "data",    Py_BuildValue("NO", PyLong_FromVoidPtr((void*)icone), Py_False),
			 "version", PyInt_FromLong(3));
  }
  %pythoncode {
  __array_interface__ = property(__array_interface__, doc='Array protocol')
  }
}
PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN
%extend Elemset {
  PyObject* __array_shape__ () {
    int m, n; self->getData(&m, &n, NULL);
    return Py_BuildValue("ii", m, n);
  }
  PyObject* __array_typestr__ () {
    char endian = PyArray_NATIVE;
    char kind   = PyArray_INTLTR;
    int  elsize = sizeof(int);
    return PyString_FromFormat("%c%c%d", endian, kind, elsize);
  }
  PyObject* __array_data__ () {
    const int* data; self->getData(NULL, NULL, &data);
    return Py_BuildValue("NO", PyString_FromFormat("%p", data), Py_False);
  }
  %pythoncode {
  __array_shape__   = property(__array_shape__,   doc='Array protocol: shape')
  __array_typestr__ = property(__array_typestr__, doc='Array protocol: typestr')
  __array_data__    = property(__array_data__,    doc='Array protocol: data')
  }
}
PYPF_NAMESPACE_END


%include "Elemset.h"

%clear  (int* nelem, int* nel);

%clear (int  n, const int elem[]);
%clear (int* n, double* elem[]);

%clear (int  nelem, int  nel, const int  icone[]);
%clear (int* nelem, int* nel, const int* icone[]);

