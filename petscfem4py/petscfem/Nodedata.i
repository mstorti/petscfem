// -*- c++ -*-
// $Id: Nodedata.i,v 1.1.2.3 2006/03/20 16:06:00 rodrigop Exp $


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

%include "Nodedata.h"

%clear (int  nnod, int* ndim, double  xyz[]);
%clear (int* nnod, int  ndim, double* xyz[]);
%clear  int* nnod, int* ndim;
