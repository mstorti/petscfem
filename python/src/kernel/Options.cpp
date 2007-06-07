// $Id$

#include "Options.h"

#include <texthash.h>

PF4PY_NAMESPACE_BEGIN

using namespace std;

Options::~Options() 
{ }

Options::Options() 
  : options(),
    hashtable(new Options::Impl)
{ }

Options::Options(Options::Impl* opt)
  : options(),
    hashtable(opt != NULL ? opt : new Options::Impl)
    
{ this->hashtable->get_entries(this->options); }

Options::Options(const Options& options) 
  : options(options.options),
    hashtable(new Options::Impl)
    

{ this->hashtable->set_entries(this->options); }


Options::Options(const map<string,string>& options) 
  : options(options),
    hashtable(new Options::Impl)
    

{ this->hashtable->set_entries(this->options); }

Options& Options::operator=(Options::Impl* options) 
{
  if (options == this->hashtable.get()) return *this;
  this->options.clear();
  this->hashtable.reset(options ? options : new Options::Impl);
  this->hashtable->get_entries(this->options);
  return *this;
}

Options& Options::operator=(const Options& opt) {
  if (&opt == this) return *this;
  this->options = opt.options;
  this->hashtable->del_entries();
  this->hashtable->set_entries(this->options);
  return *this;
}

Options& Options::operator=(const map<string,string>& options) {
  this->options = options;
  this->hashtable->del_entries();
  this->hashtable->set_entries(this->options);
  return *this;
}

bool
Options::has(const string& key) const
{
  return this->options.find(key) != this->options.end();
}

std::string
Options::get(const string& key) const
{
  map<string,string>::const_iterator m = this->options.find(key);
  if (m == this->options.end()) throw Error("option not found");
  return m->second;
}

std::string
Options::get(const string& key, const string& defval) const
{
  map<string,string>::const_iterator m = this->options.find(key);
  if (m != this->options.end()) 
    return m->second;
  else
    return defval;
}


void
Options::set(const string& key, const string& value)
{
  this->options[key] = value;
  this->hashtable->set_entry(key.c_str(), value.c_str());
}

void
Options::del(const string& key)
{
  map<string,string>::iterator m = this->options.find(key);
  if (m == this->options.end()) throw Error("option not found");
  this->options.erase(m);
  this->hashtable->del_entries();
  this->hashtable->set_entries(this->options);
}

void
Options::add(const map<string,string>& options)
{
  map<string,string>::const_iterator m = options.begin();
  while (m != options.end())
    if (not this->has(m->first)) 
      this->set(m->first, m->second);
}

void
Options::update(const map<string,string>& options)
{
  this->options.insert(options.begin(), options.end());
  this->hashtable->set_entries(options);
}

void
Options::clear()
{
  this->options.clear();
  this->hashtable->del_entries();
}

PF4PY_NAMESPACE_END



PF4PY_NAMESPACE_BEGIN

Options Options::GLOBALS;

void
Options::initGlobals() {
  Options::Impl* global_options = Options::GLOBALS;
  global_options->register_name("global_options");
  global_options->set_as_global();
  GLOBAL_OPTIONS = global_options;
}

PF4PY_NAMESPACE_END
