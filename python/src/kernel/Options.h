// $Id$

#ifndef PF4PY_OPTIONS_H
#define PF4PY_OPTIONS_H

#include <memory>
#include <string>
#include <map>
#include "petscfem4py.h"

#include "Object.h"

PF4PY_NAMESPACE_BEGIN

class Options 
  : public Object
{

public:
  typedef ::TextHashTable Impl;

protected:
  std::map<std::string,std::string> options;
  std::auto_ptr< Options::Impl >    hashtable;

#if !defined(SWIG)
public:
  typedef std::map<std::string,std::string> MapSS;
  operator Options::Impl*() const { return this->hashtable.get(); }
  operator MapSS&()             { return this->options;  }
  operator const MapSS&() const { return this->options;  }
public:
  Options(Options::Impl*);
  Options& operator=(Options::Impl*);
#endif

public:
  virtual ~Options();
  Options();
  Options(const Options& options);
  Options(const std::map<std::string,std::string>& options);
  Options& operator=(const Options&);
  Options& operator=(const std::map<std::string,std::string>&);

public:
  static Options GLOBALS;
  static void initGlobals();

public:
  bool        has(const std::string& key) const;
  std::string get(const std::string& key) const;
  std::string get(const std::string& key, const std::string& defval) const;
  void        set(const std::string& key, const std::string& value);
  void        del(const std::string& key);
  void        add(const std::map<std::string,std::string>& options);
  void        update(const std::map<std::string,std::string>& options);
  int         size()  const { return this->options.size();  }
  bool        empty() const { return this->options.empty(); }
  void        clear();

};

PF4PY_NAMESPACE_END

#endif // PF4PY_OPTIONS_H

// Local Variables:
// mode: C++
// End:
