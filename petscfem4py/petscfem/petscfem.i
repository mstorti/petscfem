// -*- c++ -*-
// $Id: petscfem.i,v 1.1.2.4 2006/03/20 16:06:00 rodrigop Exp $ 

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

%include petscinit.i
