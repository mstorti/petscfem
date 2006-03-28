// -*- c++ -*-
// $Id: petscfem.i,v 1.1.2.5 2006/03/28 22:13:25 rodrigop Exp $ 

#if 0
%module petscfem
#else
%define DOCSTRING "PESTSc-FEM for Python" %enddef
%module(package="petscfem4py",docstring=DOCSTRING) petscfem
#endif

%{
#include "petscfem4py.h"
#include "Nodedata.h"
#include "Elemset.h"
#include "Mesh.h"
#include "DofMap.h"
#include "Problem.h"
#include "NavierStokes.h"
%}

%pythoncode %{
__all__ = ['Object',
	   'Nodedata',
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
%include Nodedata.i
%include Elemset.i
%include Mesh.i
%include DofMap.i

%include Problem.i
%include NavierStokes.i
