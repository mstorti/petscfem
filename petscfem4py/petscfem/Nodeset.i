// -*- c++ -*-
// $Id: Nodeset.i,v 1.1.2.8 2006/08/22 21:53:06 dalcinl Exp $


%include Object.i

%typemap(doc, name="xnod",type="double[]") 
  (int nnod, int ndim, const double xnod[])
  "xnod: double[] value)";
ARRAY_2D(int nnod, int ndim, const double xnod[], ARRAY_INPUT, PyPF_FLOAT)
ARRAY_TYPECHECK_SEQUENCE((int nnod, int ndim, const double xnod[]), ARRAY_TYPECHECK_FLOAT)

%typemap(doc, name="array",type="double[]") 
  (int rows, int cols, const double array[])
  "array: double[] value,";
ARRAY_2D(int rows, int cols, const double array[], ARRAY_INPUT, PyPF_FLOAT)
ARRAY_TYPECHECK_SEQUENCE((int rows, int cols, const double array[]), ARRAY_TYPECHECK_FLOAT)
ARRAY_2D_NEW(int* rows, int* cols, const double* array[], PyPF_FLOAT)

%typemap(doc,name="node",type="double[]") 
  (int n, const double node[])
  "node: double[] value";
ARRAY_FLAT(int n, const double node[], ARRAY_INPUT, PyPF_FLOAT)
ARRAY_TYPECHECK_SEQUENCE((int n, const double node[]),  ARRAY_TYPECHECK_FLOAT)
ARRAY_1D_NEW(int* n, const double* node[], PyPF_FLOAT)

%apply (int*, int*, int*) { (int* nnod, int* ndim, int* nval) };

// data array accessor
PYPF_NAMESPACE_BEGIN
%extend Nodeset {
  void getArray(int *rows, int *cols, const double *array[]) {
    int  n, d, v;    self->getSizes(&n, &d, &v);
    const double *a; self->getArray(&a);
    *rows = n; *cols = d + v; *array = a;
  }
  void setArray(int rows, int cols, const double  array[]) {
    int d = self->getDim();
    self->setSizes(rows, d, cols - d);
    self->setArray(array);
  }
}
%ignore Nodeset::getArray;
%ignore Nodeset::setArray;
PYPF_NAMESPACE_END

// sequence methods
PYPF_NAMESPACE_BEGIN
%extend Nodeset {
  int  __len__()
    { return self->getSize(); }
  void __getitem__(int i, int* n, const double* node[]) 
    { self->getNode(i, n, node); }
  void __setitem__(int i, int  n, const double  node[]) 
    { self->setNode(i, n , node); }
  %pythoncode {
  def __iter__(self):
      """__iter__(self) -> iterator"""
      for i in xrange(len(self)):
          yield self[i]
  }
}
PYPF_NAMESPACE_END

// array interface #3 <http://numeric.scipy.org/array_interface.html>
PYPF_NAMESPACE_BEGIN
%extend Nodeset {
  PyObject* __array_interface__ () {
    int nnod, ndim, nval;
    const double* array; 
    self->getSizes(&nnod, &ndim, &nval);
    self->getArray(&array);
    char endian = PyArray_NATIVE;
    char kind   = PyArray_DOUBLELTR;
    int  elsize = sizeof(double);
    return Py_BuildValue("{sNsNsNsN}",
			 "shape",   Py_BuildValue("ii", nnod, ndim+nval),
			 "typestr", PyString_FromFormat("%c%c%d", endian, kind, elsize),
			 "data",    Py_BuildValue("NO", PyLong_FromVoidPtr((void*)array), Py_False),
			 "version", PyInt_FromLong(3));
  }
  %pythoncode {
  __array_interface__ = property(__array_interface__, doc='Array protocol')
  }
}
PYPF_NAMESPACE_END

// array interface #2 <http://numeric.scipy.org/array_interface.html>
PYPF_NAMESPACE_BEGIN
%extend Nodeset {
  PyObject* __array_shape__ () {
    int n, d, v; self->getSizes(&n, &d, &v);
    return Py_BuildValue("ii", n, d+v);
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
    const double* array; self->getArray(&array);
    return Py_BuildValue("NO", PyString_FromFormat("%p", array), Py_False);
  }
  %pythoncode {
  __array_shape__   = property(__array_shape__,   doc='Array protocol: shape')
  __array_typestr__ = property(__array_typestr__, doc='Array protocol: typestr')
  __array_data__    = property(__array_data__,    doc='Array protocol: data')
  }
}
PYPF_NAMESPACE_END


// properties <http://users.rcn.com/python/download/Descriptor.htm>
PYPF_NAMESPACE_BEGIN
%extend Nodeset {
  %pythoncode {
  ndim = dim   = property(getDim,  doc='number of dimensions')
  nnod = size  = property(getSize, doc='number of nodes')
  sizes = property(getSizes, doc='data layout')
  array = property(getArray, setArray, doc='data array')
  }
}
PYPF_NAMESPACE_END



%include "Nodeset.h"

%clear  (int* nnod, int* ndim, int* nval);

%clear (int  n, const double  node[]);
%clear (int* n, const double* node[]);

%clear (int  rows, int  cols, const double  array[]);
%clear (int* rows, int* cols, const double* array[]);
