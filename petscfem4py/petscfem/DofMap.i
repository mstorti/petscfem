// -*- c++ -*-
// $Id: DofMap.i,v 1.1.2.5 2006/04/27 19:09:17 rodrigop Exp $

%include Object.i
%include Mesh.i
%include Amplitude.i

PYPF_NAMESPACE_BEGIN

%pythonprepend DofMap::addFixations
%{args = len(args)==4 and (args[0],args[1:4], None) or (args[0],args[1:4], args[4])%}

%pythonprepend DofMap::addConstraints 
%{args = (args[0],args[1:4])%}

PYPF_NAMESPACE_END


ARRAY_TRIAD(int n, const int node[], const int field[], const double value[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)

ARRAY_TRIAD(int n, const int node[], const int field[], const double coeff[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)

%apply (int*, int*) { (int* local, int* global) };
%apply (int*, int*) { (int* first, int* last)   };

ARRAY_1D_NEW(int* rsize, int* ranges[], PyPF_INT)



ARRAY_FLAT(int nstt, const double stt[], ARRAY_INPUT, PyPF_FLOAT)
ARRAY_FLAT(int nsol, const double sol[], ARRAY_INPUT, PyPF_FLOAT)
ARRAY_FLAT(int nstt,       double stt[], ARRAY_OUTPUT, PyPF_FLOAT)
ARRAY_FLAT(int nsol,       double sol[], ARRAY_OUTPUT, PyPF_FLOAT)

%typemap(doc,name="state",type="double array[]") 
  (int nstt, double stt[])
  "state: double array[] value (output)";
%typemap(doc,name="state",type="double array[]") 
  (int nstt, const double stt[])
  "state: double array[] value (input)";
%typemap(doc,name="solution",type="double array[]") 
  (int nsol, double sol[])
  "solution: double array[] value (output)";
%typemap(doc,name="solution",type="double array[]") 
  (int nsol, const double sol[])
  "solution: double array[] value (input)";


%include "DofMap.h"


%clear (int n, const int node[], const int field[], const double value[]);
%clear (int n, const int node[], const int field[], const double coeff[]);

%clear (int* local, int* global);
%clear (int* first, int* last);
%clear (int* rsize, int* ranges[]);

%clear (int nstt, const double stt[]);
%clear (int nsol, const double sol[]);
%clear (int nstt,       double stt[]);
%clear (int nsol,       double sol[]);
