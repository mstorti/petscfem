// -*- c++ -*-
// $Id: Elemset.i,v 1.1.2.2 2006/03/02 21:37:12 rodrigop Exp $

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

%apply int* OUTPUT {int* nelem, int* nel};

%include "Elemset.h"
