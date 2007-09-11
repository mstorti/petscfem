// -*- c++ -*-
// $Id$

PF4PY_NAMESPACE_BEGIN
%newobject Domain::getOptions;
%newobject Domain::getNodedata;
%newobject Domain::getField;
%newobject Domain::getElemset;
%newobject Domain::getMesh;
%newobject Domain::getDofset;
%newobject Domain::getAppCtx;
PF4PY_NAMESPACE_END


%typemap(doc,name="nodes", type="int const []")    
(int nn, const int    nodes[])  "$2_name: $2_type value";
%typemap(doc,name="nodes1",type="int const []")    
(int n1, const int    nodes1[]) "$2_name: $2_type value";
%typemap(doc,name="nodes2",type="int const []")    
(int n2, const int    nodes2[]) "$2_name: $2_type value";
%typemap(doc,name="fields",type="int const []")    
(int nf, const int    fields[]) "$2_name: $2_type value";
%typemap(doc,name="values",type="double const []") 
(int nv, const double values[]) "$2_name: $2_type value";
%typemap(doc,name="coeffs",type="double const []") 
(int nc, const double coeffs[]) "$2_name: $2_type value";
ARRAY_FLAT(int nn, const int    nodes[],  ARRAY_INPUT, PyPF_INT)
ARRAY_FLAT(int n1, const int    nodes1[], ARRAY_INPUT, PyPF_INT)
ARRAY_FLAT(int n2, const int    nodes2[], ARRAY_INPUT, PyPF_INT)
ARRAY_FLAT(int nf, const int    fields[], ARRAY_INPUT, PyPF_INT)
ARRAY_FLAT(int nv, const double values[], ARRAY_INPUT, PyPF_FLOAT)
ARRAY_FLAT(int nc, const double coeffs[], ARRAY_INPUT, PyPF_FLOAT)

%include property.i
//%property(PETScFEM::Domain, PETScFEM::Comm&, comm, "communicator", self->getComm() );
%property(PETScFEM::Domain, int, ndim, "ndim", self->getNDim() );
%property(PETScFEM::Domain, int, nnod, "nnod", self->getNNod() );
%property(PETScFEM::Domain, int, ndof, "ndof", self->getNDof() );

PETSC_OBJECT_TYPEMAP(Vec)
PETSC_OBJECT_TYPEMAP(Mat)

//%apply Vec OPTIONAL { Vec x1 };
//%apply Vec OPTIONAL { Vec x0 };
%apply Vec OPTIONAL { Vec r  };
%apply Mat OPTIONAL { Mat J  };

%include "Domain.h"

%clear (int nn, const int    nodes[]);
%clear (int n1, const int    nodes1[]);
%clear (int n2, const int    nodes2[]);
%clear (int nf, const int    fields[]);
%clear (int nv, const double values[]);
%clear (int nc, const double coeffs[]);

//%clear Vec x1, Vec x0;
%clear Vec r, Mat J;
