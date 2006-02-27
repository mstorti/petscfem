// -*- c++ -*-

%{
#include "Elemset.h"
%}

ARRAY_2D(int nelem, int nel, int icone[],
	 ARRAY_INPUT, PyPF_INT)
                                                           
%typemap(check) (int nelem, int nel, int icone[])
{
  ARRAY_check_ndim(array$argnum, 2)
    //ARRAY_check_dim(array$argnum, 1, 3)
}

ARRAY_TYPECHECK_SEQUENCE((int nelem, int nel, int icone[]),
			 ARRAY_TYPECHECK_INT32)


ARRAY_2D_NEW(int* nelem, int* nel, int* icone[], PyPF_INT)

//%apply int* OUTPUT {int* size, int* nel};

%include "Elemset.h"
