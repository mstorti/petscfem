//__INSERT_LICENSE__
// $Id: autostr.cpp,v 1.9 2003/07/02 23:22:19 mstorti Exp $
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
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
int AutoString::is_conformed() const {
  for (int j=n; j>=0; j--) {
    if (!s[j]) return 1;
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void AutoString::conform() {
  if (!is_conformed()) 
    resize(n+1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
AutoString & AutoString::print() {
  printf("%s",str());
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
AutoString & AutoString::deblank(void) {
  // String should have at least the null terminator
  assert(n>0);
  // skip trailing white-space 
  int j1=0;
  while (j1<n && s[j1]!='\0' && 
         (s[j1]==' ' || s[j1]=='\t')) j1++;
  if (j1==n-1) {
    // String is completely blank, simply clear it
    assert(s[j1]=='\0');
    clear();
  } else if (s[j1]=='\"') {
    // If first non-blank character is " then
    // we take all the string from the following
    // to the previous of the matching "
    j1++;
    int j2=j1;
    while (j2<n && s[j2]!='\0' && s[j2]!='"') j2++;
    if (j2>=n || s[j2]!='"' || j2<=j1) {
      // String is empty
      clear();
    } else {
      int m=j2-j1;
      // Copy last part of string at beginning
      if (j1>0) 
      for (int k=0; k<m; k++) s[k] = s[j1+k];
      s[m]='\0';
    }
  } else {
    // Search for end of string
    int j2=j1+1;
    while (s[j2]!='\0' && s[j2]!=' ' && s[j2]!='\t') j2++;
    int m=j2-j1;
    // Copy last part of string at beginning
    if (j1>0) 
      for (int k=0; k<m; k++) s[k] = s[j1+k];
    s[m]='\0';
  }
  return *this;
}
