// -*- c++ -*-
// $Id: tpmaps.i,v 1.1.2.1 2006/04/27 19:09:18 rodrigop Exp $


// pair of int values
%typemap(in, numinputs=0, noblock=1) (int*, int*)
($*1_ltype temp1=0, $*2_ltype temp2=0)
{ $1= &temp1; $2= &temp2; }
%typemap(argout, noblock=1) (int*, int*)
{ %append_output(Py_BuildValue("ii", *$1, *$2)) ;}


// PETSc objects
%define PETSC_OBJECT_TYPEMAP(Type)
typedef struct _p_##Type* Type;
%types(Type);
%typemap(typecheck) Type = SWIGTYPE*;
%typemap(in)        Type = SWIGTYPE*;
%typemap(arginit)   Type "$1 = PETSC_NULL;";
%typemap(check, noblock=1) Type 
{ if ($1 == PETSC_NULL) %argument_nullref($type, $symname, $argnum); }
%enddef
