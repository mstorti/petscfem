// -*- c++ -*-
// $Id: Domain.i,v 1.1.2.7 2006/07/26 23:32:22 dalcinl Exp $ 

%include Object.i
%include Mesh.i
%include DofMap.i

%apply (int*, int*) { (int* nnod,  int* ndof)   };
%apply (int*, int*) { (int* local, int* global) };
%apply (int*, int*) { (int* first, int* last)   };
ARRAY_1D_NEW(int* rsize, int* ranges[], PyPF_INT)

%apply (int*, int*) { (int* start, int* end) };
ARRAY_STDVEC_OUTPUT(std::vector<int>& gdofs, PyPF_INT)
ARRAY_STDVEC_OUTPUT(std::vector<int>& ldofs, PyPF_INT)

%typemap(in,numinputs=0, noblock=1) 
std::vector<std::vector<int> >& split ($*ltype temp) {
  $1 = &temp;
}
%typemap(argout, noblock=1) 
std::vector<std::vector<int> >& split {
  for (int i=0; i<$1->size(); i++) {
    std::vector<int>& dofs = (*$1)[i];
    PyObject* o = ARRAY_NEW(&dofs[0], PyPF_INT, 1, dofs.size());
    %append_output(o);
  }
}

PYPF_NAMESPACE_BEGIN
%newobject Domain::getNodeset;
%newobject Domain::getDofset;
%newobject Domain::getMesh;
%newobject Domain::getDofMap;
PYPF_NAMESPACE_END

// properties
PYPF_NAMESPACE_BEGIN
%extend Domain {
  %pythoncode {
  nodeset = property(getNodeset, doc='Nodeset instance')
  dofset  = property(getDofset,  doc='Dofset instance')
  mesh    = property(getMesh,    doc='Mesh instance')
  dofmap  = property(getDofMap,  doc='DofMap instance')
  }
}
PYPF_NAMESPACE_END

%include "Domain.h"

%clear (int* nnod,  int* ndof);
%clear (int* local, int* global);
%clear (int* first, int* last);
%clear (int* rsize, int* ranges[]);

%clear (int* start, int* end);
%clear std::vector<int>& gdofs;
%clear std::vector<int>& ldofs;

%clear std::vector<std::vector<int> >& split;
