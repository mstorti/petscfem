// -*- c++ -*-
// $Id: Domain.i,v 1.1.2.3 2006/05/30 20:05:49 dalcinl Exp $ 

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


PYPF_NAMESPACE_BEGIN
%newobject Domain::getMesh;
%newobject Domain::getDofMap;
PYPF_NAMESPACE_END


%include "Domain.h"

%clear (int* nnod,  int* ndof);
%clear (int* local, int* global);
%clear (int* first, int* last);
%clear (int* rsize, int* ranges[]);

%clear (int* start, int* end);
%clear std::vector<int>& gdofs;
%clear std::vector<int>& ldofs;
