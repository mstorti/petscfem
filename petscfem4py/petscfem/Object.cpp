// $Id: Object.cpp,v 1.1.2.1 2006/03/06 16:56:04 rodrigop Exp $

#include "Object.h"
#include <texthash.h>

bool
PyPF::Object::hasOption(const std::string& key)
{
  OptionTable* options = this->get_opt_table();
  if (options == NULL) throw Error("null option table");
  const char* value = NULL;
  options->get_entry(key.c_str(), value);
  return (value == NULL) ? false : true;
}

std::string 
PyPF::Object::getOption(const std::string& key) 
{
  OptionTable* options = this->get_opt_table();
  if (options == NULL) throw Error("null option table");
  const char* value = NULL;
  options->get_entry(key.c_str(), value);
  if (value == NULL) throw Error("option not found");
  return value;
}

void
PyPF::Object::setOption(const std::string& key,
			const std::string& value)
{
  OptionTable* options = this->get_opt_table(true);
  if (options == NULL) throw Error("null option table");
  options->set_entry(key.c_str(), value.c_str());
}

void        
PyPF::Object::setOptions(const std::map<std::string,std::string>& M)
{
  OptionTable* options = this->get_opt_table(true);
  if (options == NULL) throw Error("null option table");
  std::map<std::string,std::string>::const_iterator m = M.begin();
  while (m != M.end()) {
    const char* key = m->first.c_str();
    const char* val = m->second.c_str(); 
    options->add_entry(key, val);
    m++;
  }
}

std::map<std::string,std::string>
PyPF::Object::getOptions()
{
  OptionTable* options = this->get_opt_table(true);
  if (options == NULL) throw Error("null option table");
  std::map<std::string,std::string> M;
  options->get_entries(M);
  return M;
}
