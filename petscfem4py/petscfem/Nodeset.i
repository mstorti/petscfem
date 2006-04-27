// -*- c++ -*-
// $Id: Nodeset.i,v 1.1.2.2 2006/04/27 19:09:17 rodrigop Exp $


%include Object.i

ARRAY_FLAT(int n, const double data[],
	   ARRAY_INPUT, PyPF_FLOAT)
ARRAY_TYPECHECK_SEQUENCE((int n, const double data[]), 
			 ARRAY_TYPECHECK_FLOAT)
ARRAY_1D_NEW(int* n, const double* data[], PyPF_FLOAT)


ARRAY_2D(int nnod, int nval, const double data[],
	 ARRAY_INPUT, PyPF_FLOAT)
ARRAY_TYPECHECK_SEQUENCE((int nnod, int nval, const double data[]), 
			 ARRAY_TYPECHECK_FLOAT)
ARRAY_2D_NEW(int* nnod, int* nval, const double* data[], PyPF_FLOAT)


%apply (int*, int*) { (int* nnod, int* nval) };


PYPF_NAMESPACE_BEGIN
%extend Nodeset {
  int  __len__() 
    { return self->getSize(); }
  void __getitem__(int i, int* n, const double* data[]) 
    { return self->getNode(i, n, data); }
  void __setitem__(int i, int  n, const double  data[]) 
    { self->setNode(i, n , data); }
  %pythoncode {
  def __iter__(self):
      """__iter__(self) -> iterator"""
      for i in xrange(len(self)):
          yield self[i]
  }
}
PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN
%extend Nodeset {
  %pythoncode {
  dim  = property(getDim,  setDim,  doc='number of dimensions')
  size = property(getSize,          doc='number of nodes')
  data = property(getData, setData, doc='nodeset data')
  }
}
PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN
%extend Nodeset {
  PyObject* __array_shape__ () {
    int m, n; self->getData(&m, &n, NULL);
    return Py_BuildValue("ii", m, n);
  }
  PyObject* __array_typestr__ () {
    PyArray_Descr *descr = PyArray_DescrFromType(PyArray_DOUBLE);
    char kind      = descr->kind;
    char byteorder = descr->byteorder;
    int  elsize    = descr->elsize;
    Py_DECREF(descr);
    return PyString_FromFormat("%c%c%d", byteorder, kind, elsize);
  }
  PyObject* __array_data__ () {
    const double* data; self->getData(NULL, NULL, &data);
    return Py_BuildValue("NO", PyString_FromFormat("%p", data), Py_False);
  }
  %pythoncode {
  __array_shape__   = property(__array_shape__,   doc='Array protocol: shape')
  __array_typestr__ = property(__array_typestr__, doc='Array protocol: typestr')
  __array_data__    = property(__array_data__,    doc='Array protocol: data')
  }
}
PYPF_NAMESPACE_END


%include "Nodeset.h"

%clear  (int* nnod, int* nval);

%clear (int  n, const double data[]);
%clear (int* n, const double* data[]);

%clear (int  nnod, int  nval, const double  data[]);
%clear (int* nnod, int* nval, const double* data[]);

