// $Id: NavierStokes.cpp,v 1.1.2.1 2006/03/02 21:37:12 rodrigop Exp $

#include "NavierStokes.h"

#include <fem.h>
#include <readmesh.h>
#include <getprop.h>
#include <utils.h>
#include <util2.h>
#include <sttfilter.h>
#include <pfmat.h>
#include <hook.h>

PyPF::NavierStokes::~NavierStokes() { }

PyPF::NavierStokes::NavierStokes(int nnod, int ndim, int ndof)
  : Problem(nnod, ndim, ndof)
{
}

void 
PyPF::NavierStokes::computeResidual() 
{ 
}

void 
PyPF::NavierStokes::computeJacobian() 
{ 
}
