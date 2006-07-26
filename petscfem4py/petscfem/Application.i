// -*- c++ -*-
// $Id: Application.i,v 1.1.2.7 2006/07/26 23:31:51 dalcinl Exp $

%include Domain.i


PYPF_NAMESPACE_BEGIN
%newobject Application::getDomain;
PYPF_NAMESPACE_END

PETSC_OBJECT_TYPEMAP(Vec)
%typemap(check, noblock=1) Vec solution "";

%include "Application.h"


PETSC_OBJECT_TYPEMAP(Vec)
%typemap(check, noblock=1) Vec x0 "";
%typemap(check, noblock=1) Vec x1 "";
%typemap(check, noblock=1) Vec r  "";

PETSC_OBJECT_TYPEMAP(Mat)
%typemap(check, noblock=1) Mat J  "";

ARRAY_STDVEC_OUTPUT(std::vector<int>& dnz, PyPF_INT)
ARRAY_STDVEC_OUTPUT(std::vector<int>& onz, PyPF_INT)

ARRAY_STDVEC_OUTPUT(std::vector<int>& xadj,   PyPF_INT)
ARRAY_STDVEC_OUTPUT(std::vector<int>& adjncy, PyPF_INT)

%include "NvrStks.h"
