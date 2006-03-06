// $Id: NavierStokes.cpp,v 1.1.2.2 2006/03/06 16:56:04 rodrigop Exp $

#include "NavierStokes.h"

#include <fem.h>
#include <readmesh.h>
#include <getprop.h>
#include <utils.h>
#include <util2.h>
#include <sttfilter.h>
#include <pfmat.h>
#include <hook.h>

PyPF::NavierStokes::~NavierStokes() 
{ }

PyPF::NavierStokes::NavierStokes() 
{ }

PyPF::NavierStokes::NavierStokes(int nnod, int ndim)
  : Problem(nnod, ndim, ndim+1)
{ }

void 
PyPF::NavierStokes::computeResidual()
{ 
}

void 
PyPF::NavierStokes::computeJacobian() 
{ 
}
