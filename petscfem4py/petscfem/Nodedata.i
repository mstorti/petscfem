// -*- c++ -*-

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

%include "Nodedata.h"
