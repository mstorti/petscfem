//__INSERT_LICENSE__
// $Id: autostr.cpp,v 1.5 2003/05/04 16:20:18 mstorti Exp $
#define _GNU_SOURCE
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <src/autostr.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString::AutoString() : s(NULL), n(0) { clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString::~AutoString() { clear(); delete[] s; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
const char *AutoString::str() const { return s; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int AutoString::size() const { return n; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::resize(int m) { 
  if (m>n) {
    char *new_s = (char *)malloc(m);
    new_s[0] = '\0';
    if (s) {
      strcpy(new_s,s);
      free(s);
    }
    s = new_s;
    n = m;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::clear() {
  resize(1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::cat(const AutoString &as) { 
  int my_len = len();
  resize(my_len + as.len()+1);
  strcpy(s+my_len, as.str());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::vsprintf(const char *f,va_list ap) { 
  n = vasprintf(&s,f,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::sprintf(const char *f,...) {
  va_list ap;
  va_start(ap,f);
  vsprintf(f,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int AutoString::len() const { return strlen(s); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::vcat_sprintf(const char *f,va_list ap) {
  AutoString t;
  t.vsprintf(f,ap);
  cat(t);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::cat_sprintf(const char *f,...) {
  va_list ap;
  va_start(ap,f);
  vcat_sprintf(f,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int AutoString::getline(FILE *fid) { 
  size_t nn = n;
  int r = ::getline(&s,&nn,fid);
  n = nn;
  return r;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::cat(const char *ss) { 
  int my_len = len();
  resize(my_len + strlen(ss)+1);
  strcpy(s+my_len,ss);
}
