// -*- c++ -*-
// $Id: Elemset.i,v 1.1.2.3 2006/03/20 16:06:00 rodrigop Exp $


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

%include "Elemset.h"

%clear (int  nelem, int  nel, int  icone[]);
%clear (int* nelem, int* nel, int* icone[]);
%clear  int* nelem, int* nel;
