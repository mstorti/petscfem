// -*- c++ -*-
// $Id: petscfem.i,v 1.1.2.11 2006/06/22 22:35:42 dalcinl Exp $ 

#if 0 // hack for numpy.distutils.command.build_ext
%module petscfem
#else
%module(package="petscfem4py", directors=1) petscfem
#endif

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

%include stl.i
%include array.i
%include macros.i
%include tpmaps.i
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
