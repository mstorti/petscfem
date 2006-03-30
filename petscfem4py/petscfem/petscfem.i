// -*- c++ -*-
// $Id: petscfem.i,v 1.1.2.6 2006/03/30 15:40:05 rodrigop Exp $ 

#if 0
%module petscfem
#else
%define DOCSTRING "PESTSc-FEM for Python" %enddef
%module(package="petscfem4py",docstring=DOCSTRING) petscfem
#endif

%{
#include "petscfem4py.h"
#include "Nodeset.h"
#include "Elemset.h"
#include "Mesh.h"
#include "DofMap.h"
#include "Problem.h"
#include "NavierStokes.h"
%}

%pythoncode %{
__all__ = ['Object',
	   'Nodeset',
	   'Elemset',
	   'Mesh',
	   'DofMap',
	   'Problem',
	   'NavierStokes']
%}

%include macros.i
%include array.i
%include stl.i

%include init.i

%import  namespace.h
%include Error.i

%include Base.i
%include Object.i
%include Nodeset.i
%include Elemset.i
%include Mesh.i
%include DofMap.i

%include Problem.i
%include NavierStokes.i
