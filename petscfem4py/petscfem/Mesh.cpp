#include "Nodedata.h"
#include "Elemset.h"
#include "Mesh.h"

#include <fem.h>

PyPF::Mesh::Mesh() 
  : Ptr(new ::Mesh) 
{ 
  (*this)->nodedata = NULL;
  (*this)->elemsetlist = da_create(sizeof(::Elemset *));
  (*this)->global_options = new TextHashTable;
}

PyPF::Mesh::~Mesh()
{ 
  /* delete[] this->ptr; */
}


// std::map<std::string,std::string> 
// PyPF::Mesh::getOptions() {
//   std::map<std::string,std::string> options;
//   return options;
// }

// void 
// PyPF::Mesh::setOptions(const std::map<std::string,std::string>& options) 
// {
//   typedef std::map<std::string,std::string> mapstr;
//   mapstr::const_iterator m = options.begin();
//   while (m!=options.end()) {
//     printf("%s : %s\n", m->first.c_str(), m->second.c_str());
//   }
// }

std::string
PyPF::Mesh::getOption(const std::string& key)
{
  const char* value = NULL;
  (*this)->global_options->get_entry(key.c_str(), value);
  return value;
}

void 
PyPF::Mesh::setOption(const std::string& key,
		      const std::string& value)
{
  (*this)->global_options->add_entry(key.c_str(), value.c_str());
}

PyPF::Nodedata
PyPF::Mesh::getNodedata()
{
  return (*this)->nodedata;
}

void 
PyPF::Mesh::setNodedata(PyPF::Nodedata& nodedata)
{
  (*this)->nodedata = nodedata;
}

void 
PyPF::Mesh::addElemset(PyPF::Elemset& elemset)
{
  ::Elemset* e = elemset;
  da_append((*this)->elemsetlist, &e);
  
}

PyPF::Elemset
PyPF::Mesh::getElemset(int i)
{
  int n = da_length((*this)->elemsetlist);
  if (i<0 || i>=n) return NULL;
  ::Elemset* e = *(::Elemset **) da_ref((*this)->elemsetlist, i);
  return e;
}

int
PyPF::Mesh::getSize()
{
  return da_length((*this)->elemsetlist);
}
  
bool
PyPF::Mesh::hasElemset(const std::string& name)
{
  ::Elemset* e = (*this)->find(name);
  return e?true:false;
}

PyPF::Elemset
PyPF::Mesh::findElemset(const std::string& name)
{
  ::Elemset* e = (*this)->find(name);
  return e;
}
