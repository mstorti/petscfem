// -*- c++ -*-
// $Id: Nodeset.i,v 1.1.2.1 2006/03/30 15:40:05 rodrigop Exp $


%include Object.i


ARRAY_2D(int nnod, int ndim, double xyz[],
	 ARRAY_INPUT, PyPF_FLOAT)
                                                           
%typemap(check) (int nnod, int ndim, double xyz[])
{
  ARRAY_check_ndim(array$argnum, 2)
}

ARRAY_TYPECHECK_SEQUENCE((int nnod, int ndim, double xyz[]), 
			 ARRAY_TYPECHECK_FLOAT)


ARRAY_2D_NEW(int* nnod, int* ndim, double* xyz[], PyPF_FLOAT)

%apply int* OUTPUT {int* nnod, int* ndim};



PYPF_NAMESPACE_BEGIN

%template() ::std::vector<double>;

%extend Nodeset {
  int __len__() { 
    int nnod;  
    self->getSize(&nnod, NULL);
    return nnod; 
  }
  Node __getitem__(int n) { 
    return self->getNode(n); 
  }
  void __setitem__(int n, const Node& node) {
    self->setNode(n, node); 
  }
}

%feature("shadow") Nodeset::__iter__ %{
def __iter__(self):
    for i in xrange(len(self)):
        yield self[i]
%}
%extend Nodeset { void __iter__() { } }

PYPF_NAMESPACE_END

%include "Nodeset.h"

%clear (int  nnod, int* ndim, double  xyz[]);
%clear (int* nnod, int  ndim, double* xyz[]);
%clear  int* nnod, int* ndim;
