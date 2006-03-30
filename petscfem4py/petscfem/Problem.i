// -*- c++ -*-
// $Id: Problem.i,v 1.1.2.3 2006/03/30 15:18:14 rodrigop Exp $

%include petsc.i

PETSC_OBJECT_TYPEMAP(Vec)

PYPF_NAMESPACE_BEGIN

%newobject Problem::getMesh;
%newobject Problem::getDofMap;

PYPF_NAMESPACE_END

%apply int* OUTPUT {int* local, int* global};
%apply int* OUTPUT {int* first, int* last};

%include "Problem.h"

%clear int* local, int* global;
%clear int* first, int* last;
