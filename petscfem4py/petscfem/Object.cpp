// $Id: Object.cpp,v 1.1.2.3 2006/03/28 22:13:25 rodrigop Exp $

#include "Object.h"

#include <petsc.h>
#include <texthash.h>


PYPF_NAMESPACE_BEGIN

using namespace std;

Object::~Object() 
{ }

Object::Object()
  : refcnt(0), comm(PETSC_COMM_WORLD), options()
{ }

Object::Object(const Object& obj)
  : refcnt(0), comm(obj.comm), options(obj.options)
{ }


MPI_Comm&
Object::getComm() const 
{
  return const_cast<MPI_Comm&>(this->comm);
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


PYPF_NAMESPACE_END
