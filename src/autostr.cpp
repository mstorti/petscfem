//__INSERT_LICENSE__
// $Id: autostr.cpp,v 1.8 2003/05/04 16:43:50 mstorti Exp $
#define _GNU_SOURCE
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <src/autostr.h>

int AutoString::total=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString::AutoString() : s(NULL), n(0) { clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString::~AutoString() { clear(); total -= n; free(s); s=NULL; n=0; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
const char *AutoString::str() const { return s; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int AutoString::size() const { return n; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::resize(int m) { 
  // Allocate new size and copy old vector
  if (m>n) {
    char *new_s = (char *)malloc(m);
    total += m;
    new_s[0] = '\0';
    if (s) {
      strcpy(new_s,s);
      free(s);
      total -= n;
    }
    s = new_s;
    n = m;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString & AutoString::clear() {
  total -= n;
  free(s);
  s = NULL;
  n = 0;
  resize(1);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString & AutoString::cat(const AutoString &as) { 
  int my_len = len();
  resize(my_len + as.len()+1);
  strcpy(s+my_len, as.str());
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString &  AutoString::vsprintf(const char *f,va_list ap) { 
  if (s) {
    free(s);
    n=0;
    s=NULL;
  }
  n = vasprintf(&s,f,ap);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString &  AutoString::sprintf(const char *f,...) {
  va_list ap;
  va_start(ap,f);
  return vsprintf(f,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int AutoString::len() const { return strlen(s); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString & AutoString::vcat_sprintf(const char *f,va_list ap) {
  AutoString t;
  t.vsprintf(f,ap);
  return cat(t);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString &  AutoString::cat_sprintf(const char *f,...) {
  va_list ap;
  va_start(ap,f);
  return vcat_sprintf(f,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int AutoString::getline(FILE *fid) { 
  size_t nn = n;
  int r = ::getline(&s,&nn,fid);
  n = nn;
  return r;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString & AutoString::cat(const char *ss) { 
  int my_len = len();
  resize(my_len + strlen(ss)+1);
  strcpy(s+my_len,ss);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString &  AutoString::set(const char *ss) { 
  clear();
  return cat(ss);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString & AutoString::set(const AutoString &as) { 
  clear();
  return cat(as);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString & AutoString::vfprintf(FILE *fid,va_list ap) {
  ::vfprintf(fid,str(),ap);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString & AutoString::fprintf(FILE *fid,...) {
  va_list ap;
  va_start(ap,fid);
  return vfprintf(fid,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
const AutoString & AutoString::print() const {
  printf("%s",str());
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString & AutoString::print() {
  printf("%s",str());
  return *this;
}
