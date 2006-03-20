// -*- c++ -*-
// $Id: Problem.i,v 1.1.2.2 2006/03/20 16:06:00 rodrigop Exp $

%include petsc.i

PETSC_OBJECT_TYPEMAP(Vec)

PYPF_NAMESPACE_BEGIN

%newobject Problem::getMesh;
%newobject Problem::getDofMap;

PYPF_NAMESPACE_END


%include "Problem.h"
