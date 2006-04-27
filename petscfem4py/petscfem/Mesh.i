// -*- c++ -*-
// $Id: Mesh.i,v 1.1.2.7 2006/04/27 19:09:17 rodrigop Exp $


%include Object.i
%include Nodeset.i
%include Elemset.i


PYPF_NAMESPACE_BEGIN
%newobject Mesh::getNodeset;
%newobject Mesh::getElemset;
PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN
%extend Mesh {
  int __len__() 
    { return self->getSize(); }
  %newobject __getitem__;
  Elemset& __getitem__(int i)
    { return self->getElemset(i); }
  void __setitem__(int i, Elemset& elemset)
    { self->setElemset(i, elemset); }
  void __delitem__(int i)
    { self->delElemset(i); }
  %pythoncode {
  def __iter__(self):
    """__iter__(self) -> iterator"""
    for i in xrange(len(self)):
        yield self[i]
  }
}
PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN
%extend Mesh {
  %pythoncode {
  size = property(getSize, doc='mesh size (number of elemsets)')
  }
}
PYPF_NAMESPACE_END


%include "Mesh.h"
