// -*- c++ -*-
// $Id: Problem.i,v 1.1.2.4 2006/04/27 19:09:17 rodrigop Exp $

%include Object.i
%include Mesh.i
%include DofMap.i

ARRAY_1D_NEW(int* n, int* dofs[], PyPF_INT)
ARRAY_1D_FREEARG(int* n, int* dofs[], delete[])


PYPF_NAMESPACE_BEGIN
%newobject Problem::getMesh;
%newobject Problem::getDofMap;
PYPF_NAMESPACE_END

%apply (int*, int*) { (int* nnod,  int* ndof) };
%apply (int*, int*) { (int* local, int* global) };
%apply (int*, int*) { (int* first, int* last)   };

PETSC_OBJECT_TYPEMAP(Vec);

%include "Problem.h"

%clear (int* nnod,  int* ndof);
%clear (int* local, int* global);
%clear (int* first, int* last);
