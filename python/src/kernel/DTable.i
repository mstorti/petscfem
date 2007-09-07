// -*- c++ -*-
// $Id$

%template() std::pair<int,int>;
%template() std::pair<char,int>;

%typemap(doc, name="array",type="int[]") 
(int rows, int cols, const int array[])
"$3_name: $3_type value";;
ARRAY_2D(int rows, int cols, const int array[], ARRAY_INPUT, PyPF_INT)
ARRAY_TYPECHECK_SEQUENCE((int rows, int cols, const int array[]), ARRAY_TYPECHECK_INT32)

%typemap(doc, name="array",type="double[]") 
(int rows, int cols, const double array[])
"$3_name: $3_type value";;
ARRAY_2D(int rows, int cols, const double array[], ARRAY_INPUT, PyPF_FLOAT)
ARRAY_TYPECHECK_SEQUENCE((int rows, int cols, const double array[]), ARRAY_TYPECHECK_FLOAT)

%include "DTable.h"

// PF4PY_NAMESPACE_BEGIN
// %feature("implicitconv") DTable<int>;
// PF4PY_NAMESPACE_END

// PF4PY_NAMESPACE_BEGIN
// %feature("implicitconv") DTable<double>;
// PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
%ignore DTable<int>::operator int*;
%ignore DTable<int>::operator const int*;
PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
%ignore DTable<double>::operator double*;
%ignore DTable<double>::operator const double*;
PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
%template(DTableI) DTable<int>;
%template(DTableS) DTable<double>;
PF4PY_NAMESPACE_END

%array_interface( PF4PY_NAMESPACE::DTable<int> );
%array_interface( PF4PY_NAMESPACE::DTable<double> );
