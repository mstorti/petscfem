// -*- c++ -*-
%module petscfem

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
__all__  = ['Int', 'Float']
__all__ += ['Nodedata',
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
%include SmartPtr.i

%include Nodedata.i
%include Elemset.i
%include Mesh.i
%include DofMap.i

%include Problem.i
%include NavierStokes.i

%include petscinit.i
