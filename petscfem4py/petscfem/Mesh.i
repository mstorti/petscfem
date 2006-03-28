// -*- c++ -*-
// $Id: Mesh.i,v 1.1.2.4 2006/03/28 22:13:25 rodrigop Exp $


%include Object.i
%include Nodedata.i
%include Elemset.i


PYPF_NAMESPACE_BEGIN

%newobject Mesh::getNodedata;
%newobject Mesh::getElemset;

%template() ::std::vector<Elemset*>;

%extend Mesh {

  int __len__() { 
    return self->getSize();
  }

  %newobject __getitem__;
  Elemset* __getitem__(int i) { 
    return self->getElemset(i);
  }

  void __setitem__(int i, Elemset* e) { 
    self->setElemset(i, e);
  }

}

%feature("shadow") Mesh::__iter__ %{
def __iter__(self):
    for i in xrange(len(self)):
        yield self[i]
%}
%extend Mesh { void __iter__() { } }

PYPF_NAMESPACE_END


%include "Mesh.h"
