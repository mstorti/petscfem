// -*- c++ -*-
// $Id: Elemset.i,v 1.1.2.5 2006/03/30 15:18:14 rodrigop Exp $


%include Object.i


ARRAY_2D(int nelem, int nel, int icone[],
	 ARRAY_INPUT, PyPF_INT)
                                                           
%typemap(check) (int nelem, int nel, int icone[])
{
  ARRAY_check_ndim(array$argnum, 2)
}

ARRAY_TYPECHECK_SEQUENCE((int nelem, int nel, int icone[]),
			 ARRAY_TYPECHECK_INT32)


ARRAY_2D_NEW(int* nelem, int* nel, int* icone[], PyPF_INT)

%apply int* OUTPUT {int* nelem, int* nel};


PYPF_NAMESPACE_BEGIN

%template() ::std::vector<int>;

%pythonappend  Elemset::getConnectivity %{val-=1%}

%extend Elemset {
  int __len__() { 
    int nelem;
    self->getSize(&nelem, NULL);
    return nelem;
  }
  Elem __getitem__(int n) { 
    return self->getElem(n); 
  }
  void __setitem__(int n, const Elem& elem) {
    self->setElem(n, elem); 
  }
}
%feature("shadow") Elemset::__iter__ %{
def __iter__(self):
    for i in xrange(len(self)):
        yield self[i]
%}
%extend Elemset { void __iter__() { } }

PYPF_NAMESPACE_END


%include "Elemset.h"

%clear (int  nelem, int  nel, int  icone[]);
%clear (int* nelem, int* nel, int* icone[]);
%clear  int* nelem, int* nel;
