// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: autostr.h,v 1.6 2003/05/04 16:28:15 mstorti Exp $
#ifndef PETSCFEM_AUTOSTR_H
#define PETSCFEM_AUTOSTR_H
#include <cstdarg>
#include <cstdio>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class AutoString {
private:
  char *s;
  int n;
public:
  AutoString();
  ~AutoString();
  const char *str() const;
  int size() const;
  int len() const;
  void resize(int m);
  AutoString & clear();
  AutoString & sprintf(const char *,...);
  AutoString & vsprintf(const char *,va_list ap);
  AutoString & cat_sprintf(const char *,...);
  AutoString & vcat_sprintf(const char *,va_list ap);
  AutoString & cat(const AutoString &s);
  AutoString & cat(const char *s);
  AutoString & vfprintf(FILE *fid,va_list ap);
  AutoString & fprintf(FILE *fid,...);
  const AutoString & print() const;
  AutoString & print();
  AutoString & set(const AutoString &s);
  AutoString & set(const char *s);
  int getline(FILE *fid);
};

#endif
