// -*- c++ -*-
// $Id: petscfem.i,v 1.1.2.12 2006/06/23 21:43:04 dalcinl Exp $ 

/*
%module petscfem // hack for numpy.distutils, do not remove !!
*/

%module(package="petscfem4py", directors=1) petscfem

%header %{
#include "petscfem4py.h"
#include "Options.h"
#include "Object.h"
#include "Nodeset.h"
#include "Elemset.h"
#include "Mesh.h"
#include "Amplitude.h"
#include "Dofset.h"
#include "DofMap.h"
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
%include Options.i
%include Comm.i
%include Object.i
%include Nodeset.i
%include Elemset.i
%include Mesh.i
%include Amplitude.i
%include Dofset.i
%include DofMap.i
%include Domain.i


%header %{
#include "Application.h"
#include "NvrStks.h"
%}
%include Application.i
