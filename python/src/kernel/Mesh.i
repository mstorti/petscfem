// -*- c++ -*-
// $Id$

PF4PY_NAMESPACE_BEGIN
PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
%newobject Mesh::getNodedata;
%newobject Mesh::getElemset;
%newobject Mesh::getOptions;
PF4PY_NAMESPACE_END

ARRAY_STDVEC_OUTPUT(std::vector<int>& npart, PyPF_INT)

%include "Mesh.h"

%clear std::vector<int>& npart;
