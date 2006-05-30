// -*- c++ -*-
// $Id: Application.i,v 1.1.2.2 2006/05/30 20:16:44 dalcinl Exp $

%include Domain.i


PYPF_NAMESPACE_BEGIN
%newobject Application::getDomain;
PYPF_NAMESPACE_END


ARRAY_FLAT(int ns, const double state[],
	   ARRAY_INPUT, PyPF_FLOAT)
ARRAY_FLAT(int nn, const int nodes[],
	   ARRAY_INPUT, PyPF_INT)
ARRAY_FLAT(int nf, const int fields[],
	   ARRAY_INPUT, PyPF_INT)

%header %{
#include <utility>
#include <vector>
%}

ARRAY_STDVEC_OUTPUT(std::vector<double>& values,  PyPF_FLOAT)

%typemap(in, numinputs=0) (std::pair<int,int>& shape, std::vector<double>& values)
  ($*1_ltype values, $*2_ltype shape) "$1 = &values, $2 = &shape;"
%typemap(argout) (std::pair<int,int>& shape, std::vector<double>& values)
{
  PyObject* o = ARRAY_NEW(&(*$2)[0], PyArray_DOUBLE, 2, $1->first, $1->second);
  ARRAY_arg_fail($symname, $argnum);
  %append_output(o);
}

PYPF_NAMESPACE_BEGIN

%extend Application {

  void getNodalValues(int ns, const double state[], double time,
		      std::vector<double>& values)
  {
    if (ns == 0) return;
    else if (ns != self->getDomain().getDofSize())
      throw PETScFEM::Error("invalid size for state vector");
    values.resize(self->getDomain().getSize());
    self->getNodalValues(state,time,&values[0]);
  }

  void getNodalValues(int ns, const double state[], double time,
		      int nn, const int nodes[],
		      int nf, const int fields[],
		      std::pair<int,int>& shape, std::vector<double>& values)
  {
    if (nn == 0 || nf == 0 || ns == 0) return;
    else if (ns != self->getDomain().getDofSize())
      throw PETScFEM::Error("invalid size for state vector");
    shape.first = nn; shape.second = nf;
    values.resize(nn*nf);
    self->getNodalValues(state,time,nn,nodes,nf,fields,&values[0]);
  }

}

%ignore Application::getNodalValues;

PYPF_NAMESPACE_END

%include "Application.h"

%clear (int ns, const double state[]);
%clear (int nn, const int nodes[]);
%clear (int nf, const int field[]);
%clear (std::vector<double>& values);
%clear (std::pair<int,int>& shape, std::vector<double>& values);




PETSC_OBJECT_TYPEMAP(Vec)
PETSC_OBJECT_TYPEMAP(Mat)

%typemap(check, noblock=1) Vec x0 "";
%typemap(check, noblock=1) Vec x1 "";
%typemap(check, noblock=1) Mat J  "";

%include "NvrStks.h"
