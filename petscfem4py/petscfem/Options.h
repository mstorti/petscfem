// $Id: Options.h,v 1.1.2.2 2006/08/22 22:10:43 dalcinl Exp $

#ifndef PYPF_OPTIONS_H
#define PYPF_OPTIONS_H

#include <string>
#include <map>
#include "petscfem4py.h"


PYPF_NAMESPACE_BEGIN

class Options {
  
  typedef ::Options Base;

protected:
  Base*                             texthash;
  std::map<std::string,std::string> options;

public:
  virtual ~Options();
  Options();
  Options(const Options& options);
  Options(const std::map<std::string,std::string>& options);

public:
  Options& operator=(const Options&);
  Options& operator=(const std::map<std::string,std::string>&);

public:
  static Options GLOBAL;
  static void init();

public:
  bool        has(const std::string& key) const;
  std::string get(const std::string& key) const;
  void        set(const std::string& key, const std::string& value);
  void        del(const std::string& key);

  void   add(const std::map<std::string,std::string>& options);
  void   update(const std::map<std::string,std::string>& options);
  int    size()  const { return this->options.size(); }
  bool   empty() const { return this->options.empty(); }
  void   clear();

#ifndef SWIG
public:
  Options(Options::Base*);
  Options& operator=(Options::Base*);
public:
  typedef std::map<std::string,std::string> Table;
  operator Base*&() { return this->texthash; }
  operator Table&() { return this->options;  }
  operator Base* const&() const { return this->texthash; }
  operator const Table&() const { return this->options;  }
#endif

};

PYPF_NAMESPACE_END

#endif // PYPF_OPTIONS_H

// Local Variables:
// mode: C++
// End:
