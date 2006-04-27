// $Id: Options.h,v 1.1.2.1 2006/04/27 19:09:17 rodrigop Exp $

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
  Options(Options::Base*);
  Options(const Options&);
  Options(const std::map<std::string,std::string>&);

public:
  Options& operator=(Options::Base*);
  Options& operator=(const Options&);
  Options& operator=(const std::map<std::string,std::string>&);

public:
  bool        has(const std::string& key) const;
  std::string get(const std::string& key) const;
  void        set(const std::string& key, const std::string& value);
  void        del(const std::string& key);
  void   set(const std::map<std::string,std::string>&);
  void   add(const std::map<std::string,std::string>&);

  void   clear();
  bool   empty() const { return this->options.empty(); }
  int    size()  const { return this->options.size();  }

#ifndef SWIG
public:
  typedef std::map<std::string,std::string> Table;
  operator Base*&() { return this->texthash; }
  operator Table&() { return this->options;  }
  operator Base* const&() const { return this->texthash; }
  operator const Table&() const { return this->options;  }
#endif

};


PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN

class OPTIONS
{

private:
  OPTIONS(const OPTIONS&);

public:
  static Options GLOBAL;
  static void init();

public:
  ~OPTIONS() { };
  OPTIONS() { };

public:
  static bool hasOption(const std::string& key)
  {  return GLOBAL.has(key); }
  static std::string getOption(const std::string& key) 
  { return GLOBAL.get(key); }
  static void setOption(const std::string& key, const std::string& value) 
  { GLOBAL.set(key, value); }
  static void delOption(const std::string& key)
  { GLOBAL.del(key); }

};

PYPF_NAMESPACE_END


#endif // PYPF_OPTIONS_H

// Local Variables:
// mode: C++
// End:
