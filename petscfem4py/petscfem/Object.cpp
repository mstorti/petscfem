// $Id: Object.cpp,v 1.1.2.5 2006/04/27 19:09:17 rodrigop Exp $

#include "Object.h"

#include <petsc.h>
#include <texthash.h>


PYPF_NAMESPACE_BEGIN

using namespace std;

Object::~Object() 
{ }

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

Comm
Object::getComm() const
{
  return this->comm;
}

void
Object::setComm(Comm& comm)
{
  PYPF_ASSERT((MPI_Comm)comm!=MPI_COMM_NULL, "cannot set null communicator");
  this->comm = comm;
}

void
Object::setComm(MPI_Comm comm)
{
  PYPF_ASSERT(comm!=MPI_COMM_NULL, "cannot set null communicator");
  this->comm = comm;
}

bool
Object::hasOption(const string& key) const
{
  return this->options.has(key);
}

string 
Object::getOption(const string& key) const
{
  return this->options.get(key);
}

void
Object::setOption(const string& key,
		  const string& value)
{
  return this->options.set(key, value);
}

map<string,string>
Object::getOptions() const
{
  return this->options;
}

void        
Object::setOptions(const map<string,string>& M)
{
  return this->options.set(M);
}

void        
Object::addOptions(const map<string,string>& M)
{
  return this->options.add(M);
}

void        
Object::delOptions()
{
  this->options.clear();
}


PYPF_NAMESPACE_END
