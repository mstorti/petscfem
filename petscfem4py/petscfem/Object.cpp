// $Id: Object.cpp,v 1.1.2.8 2006/08/22 22:10:43 dalcinl Exp $

#include "Object.h"

#include <petsc.h>
#include <texthash.h>


PYPF_NAMESPACE_BEGIN

using namespace std;

Object::~Object() 
{ this->comm = MPI_COMM_NULL; }

Object::Object()
  : RefCounter(),
    comm(PETSC_COMM_WORLD), options()
{ }

Object::Object(const Object& obj)
  : RefCounter(obj),
    comm(obj.comm), options(obj.options)
{ }

Object::Object(MPI_Comm comm)
  : RefCounter(),
    comm(comm), options()
{ }

MPI_Comm
Object::getComm() const
{
  return this->comm;
}

void
Object::setComm(MPI_Comm comm)
{
  PYPF_ASSERT(comm!=MPI_COMM_NULL, "cannot set null communicator");
  this->comm = comm;
}

PYPF_NAMESPACE_END
