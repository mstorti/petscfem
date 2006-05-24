// -*- c++ -*-
// $Id: Application.i,v 1.1.2.1 2006/05/24 21:02:33 dalcinl Exp $

%include Domain.i


PYPF_NAMESPACE_BEGIN
%newobject Application::getDomain;
PYPF_NAMESPACE_END
  
%include "Application.h"



PETSC_OBJECT_TYPEMAP(Vec)
PETSC_OBJECT_TYPEMAP(Mat)

%typemap(check, noblock=1) Vec x0 "";
%typemap(check, noblock=1) Vec x1 "";
%typemap(check, noblock=1) Mat J  "";

%include "NvrStks.h"
