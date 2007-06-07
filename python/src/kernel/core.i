// -*- c++ -*-
// $Id$ 

/* ---------
%module core // hack for numpy.distutils, do not remove !!
------------ */

%module(package="pf4py.kernel", directors=1) core

%feature("python:nondynamic", 0);

%header %{
#include "namespace.h"

#include "Comm.h"

#include "Object.h"

#include "Options.h"
#include "DTable.h"
#include "PTable.h"

#include "Elemset.h"
#include "Mesh.h"

#include "Amplitude.h"
#include "Dofset.h"

#include "AppCtx.h"
#include "Domain.h"
%}

%include tpmaps.i
%include stl.i
%include array.i
%include macros.i
%import  namespace.h

%include Docs.i
%include Init.i
%include Error.i

%ignore *::operator=;

%include Comm.i

%include Object.i

%include Options.i
%include DTable.i
%include PTable.i

%include Elemset.i
%include Mesh.i

%include Amplitude.i
%include Dofset.i

%include AppCtx.i
%include Domain.i
