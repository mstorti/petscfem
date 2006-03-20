// -*- c++ -*-
// $Id: Mesh.i,v 1.1.2.3 2006/03/20 16:06:00 rodrigop Exp $


%include Object.i
%include Nodedata.i
%include Elemset.i


PYPF_NAMESPACE_BEGIN

%newobject Mesh::getNodeData;
%newobject Mesh::getElemset;

%template() ::std::vector<Elemset*>;

PYPF_NAMESPACE_END


%include "Mesh.h"
