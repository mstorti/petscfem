// -*- c++ -*-
// $Id: Domain.i,v 1.1.2.1 2006/05/24 21:04:02 dalcinl Exp $ 

%include Object.i
%include Mesh.i
%include DofMap.i

ARRAY_1D_NEW(int* n, int* dofs[], PyPF_INT)
ARRAY_1D_FREEARG(int* n, int* dofs[], delete[])

%apply (int*, int*) { (int* nnod,  int* ndof) };
%apply (int*, int*) { (int* local, int* global) };
%apply (int*, int*) { (int* first, int* last)   };

ARRAY_STDVEC_OUTPUT(std::vector<int>& ldofs,  PyPF_INT)
ARRAY_STDVEC_OUTPUT(std::vector<int>& xadj,   PyPF_INT)
ARRAY_STDVEC_OUTPUT(std::vector<int>& adjncy, PyPF_INT)


PYPF_NAMESPACE_BEGIN
%newobject Domain::getMesh;
%newobject Domain::getDofMap;
PYPF_NAMESPACE_END

%include "Domain.h"

%clear (int* nnod,  int* ndof);
%clear (int* local, int* global);
%clear (int* first, int* last);
