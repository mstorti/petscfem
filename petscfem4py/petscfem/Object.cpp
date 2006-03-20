// $Id: Object.cpp,v 1.1.2.2 2006/03/20 16:06:00 rodrigop Exp $

#include "Object.h"

#include <petsc.h>
#include <texthash.h>


PYPF_NAMESPACE_BEGIN

OptionTable*
Object::get_opt_table() const
{ 
  return NULL; 
}

Object::~Object() 
{ 
  //PYPF_DELETE_SCLR(this->options);
}

Object::Object()
  : refcnt(0), comm(PETSC_COMM_WORLD), options()
{ }

Object::Object(const Object& obj)
  : refcnt(0), comm(PETSC_COMM_WORLD), options(obj.options)
{ }


MPI_Comm&
Object::getComm() const 
{
  return const_cast<MPI_Comm&>(this->comm);
}

bool
Object::hasOption(const std::string& key) const
{
  OptionTable* options = this->get_opt_table();
  if (options == NULL) throw Error("null pointer for option table");
  const char* value = NULL;
  options->get_entry(key.c_str(), value);
  return (value == NULL) ? false : true;
}

std::string 
Object::getOption(const std::string& key) const
{
  OptionTable* options = this->get_opt_table();
  if (options == NULL) throw Error("null pointer for option table");
  const char* value = NULL;
  options->get_entry(key.c_str(), value);
  if (value == NULL) throw Error("option not found");
  return value;
}

void
Object::setOption(const std::string& key,
			const std::string& value)
{
  OptionTable* options = this->get_opt_table();
  if (options == NULL) throw Error("null pointer for option table");
  options->set_entry(key.c_str(), value.c_str());
}

void        
Object::setOptions(const std::map<std::string,std::string>& M)
{
  OptionTable* options = this->get_opt_table();
  if (options == NULL) throw Error("null pointer for option table");
  std::map<std::string,std::string>::const_iterator m = M.begin();
  while (m != M.end()) {
    const char* key = m->first.c_str();
    const char* val = m->second.c_str(); 
    options->add_entry(key, val);
    m++;
  }
}

std::map<std::string,std::string>
Object::getOptions() const
{
  OptionTable* options = this->get_opt_table();
  if (options == NULL) throw Error("null pointer for option table");
  std::map<std::string,std::string> M;
  options->get_entries(M);
  return M;
}

PYPF_NAMESPACE_END
