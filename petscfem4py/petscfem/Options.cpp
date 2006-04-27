// $Id: Options.cpp,v 1.1.2.1 2006/04/27 19:09:17 rodrigop Exp $

#include "Options.h"

#include <texthash.h>


PYPF_NAMESPACE_BEGIN

using namespace std;

Options::~Options() 
{ PYPF_DELETE_SCLR(this->texthash); }

Options::Options() 
  : texthash(new Options::Base),
    options()
{ }

Options::Options(Options::Base* opt)
  : texthash(opt != NULL ? opt : new Options::Base),
    options()
{ this->texthash->get_entries(this->options); }

Options::Options(const Options& options) 
  : texthash(new Options::Base),
    options(options.options)

{ this->texthash->set_entries(this->options); }


Options::Options(const map<string,string>& options) 
  : texthash(new Options::Base),
    options(options)

{ this->texthash->set_entries(this->options); }

Options& Options::operator=(Options::Base* options) 
{
  if (options == this->texthash) return *this;
  this->options.clear();
  PYPF_DELETE_SCLR(this->texthash);
  this->texthash = options ? options : new Options::Base;
  this->texthash->get_entries(this->options);
}

Options& Options::operator=(const Options& opt) {
  if (&opt == this) return *this;
  this->options = opt.options;
  this->texthash->del_entries();
  this->texthash->set_entries(this->options);
  return *this;
}

Options& Options::operator=(const map<string,string>& options) {
  this->options = options;
  this->texthash->del_entries();
  this->texthash->set_entries(this->options);
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

void
Options::set(const string& key, const string& value)
{
  this->options[key] = value;
  this->texthash->set_entry(key.c_str(), value.c_str());
}

void
Options::del(const string& key)
{
  map<string,string>::iterator m = this->options.find(key);
  if (m == this->options.end()) throw Error("option not found");
  this->options.erase(m);
  this->texthash->del_entries();
  this->texthash->set_entries(this->options);
}

void
Options::set(const map<string,string>& options)
{
  this->options = options;
  this->texthash->del_entries();
  this->texthash->set_entries(this->options);
}

void
Options::add(const map<string,string>& options)
{
  this->options.insert(options.begin(), options.end());
  this->texthash->set_entries(options);
}

void
Options::clear()
{
  this->options.clear();
  this->texthash->del_entries();
}

PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN

Options OPTIONS::GLOBAL;

void
OPTIONS::init() {
  TextHashTable* global_options = OPTIONS::GLOBAL;
  global_options->register_name("global_options");
  global_options->set_as_global();
  GLOBAL_OPTIONS = global_options;
}
PYPF_NAMESPACE_END
