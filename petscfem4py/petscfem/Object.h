// -*- c++ -*-

#ifndef PYPF_OBJECT_H
#define PYPF_OBJECT_H

#include <string>
#include <map>
#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN

class Object
{

protected:
  virtual OptionTable* get_opt_table(bool create=false) = 0;

public:
  bool        hasOption(const std::string& key);
  std::string getOption(const std::string& key);
  void        setOption(const std::string& key,
			const std::string& value);

  std::map<std::string,std::string> getOptions();
  void setOptions(const std::map<std::string,std::string>&);
};


PYPF_NAMESPACE_END

#endif // PYPF_OBJECT_H
