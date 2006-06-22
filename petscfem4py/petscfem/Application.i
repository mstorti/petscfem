// -*- c++ -*-
// $Id: Application.i,v 1.1.2.5 2006/06/22 22:35:42 dalcinl Exp $

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

%include "NvrStks.h"
