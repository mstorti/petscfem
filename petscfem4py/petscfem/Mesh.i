// -*- c++ -*-
// $Id: Mesh.i,v 1.1.2.9 2006/06/26 19:43:19 dalcinl Exp $


%include Object.i
%include Nodeset.i
%include Elemset.i


PYPF_NAMESPACE_BEGIN
%newobject Mesh::getNodeset;
%newobject Mesh::getElemset;
PYPF_NAMESPACE_END


// special methods
PYPF_NAMESPACE_BEGIN
%extend Mesh {
  int __len__() 
    { return self->getSize(); }
  %newobject __getitem__;
  Elemset& __getitem__(int i)
    { return self->getElemset(i); }
  %pythoncode {
  def __iter__(self):
    """__iter__(self) -> iterator"""
    for i in xrange(len(self)):
        yield self[i]
  }
}
PYPF_NAMESPACE_END


// properties
PYPF_NAMESPACE_BEGIN
%extend Mesh {
  %pythoncode {
  nodeset = property(getNodeset, doc='nodeset instance')
  size    = property(getSize,    doc='number of elemset instances')
  }
}
PYPF_NAMESPACE_END


%include "Mesh.h"
