//__INSERT_LICENSE__
// $Id: autostr.cpp,v 1.2 2003/02/16 15:15:53 mstorti Exp $
#define _GNU_SOURCE
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <src/autostr.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString::AutoString() : s(NULL), n(0) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString::~AutoString() { clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
const char *AutoString::str() const { return s; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int AutoString::size() { return n; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::resize(int m) { 
  if (m>n) {
    char *new_s = (char *)malloc(m);
    if (n>0) {
      strcpy(new_s,s);
      free(s);
      s = new_s;
    }
    n=m;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::clear() {
    if (n>0) { n=0; free(s); s=NULL; }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::cat(AutoString &as) { 
  int my_len = len();
  resize(my_len + as.len()+1);
  strcpy(s+my_len, as.str());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::vsprintf(char *f,va_list ap) { 
  n = vasprintf(&s,f,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::sprintf(char *f,...) {
  va_list ap;
  va_start(ap,f);
  vsprintf(f,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int AutoString::len() { return strlen(s); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::vcat_sprintf(char *f,va_list ap) {
  AutoString t;
  t.vsprintf(f,ap);
  cat(t);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AutoString::cat_sprintf(char *f,...) {
  va_list ap;
  va_start(ap,f);
  vcat_sprintf(f,ap);
}
