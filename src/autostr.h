// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: autostr.h,v 1.5 2003/05/04 16:20:18 mstorti Exp $
#ifndef PETSCFEM_AUTOSTR_H
#define PETSCFEM_AUTOSTR_H
#include <cstdarg>

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
  void clear();
  void sprintf(const char *,...);
  void vsprintf(const char *,va_list ap);
  void cat_sprintf(const char *,...);
  void vcat_sprintf(const char *,va_list ap);
  void cat(const AutoString &s);
  void cat(const char *s);
  int getline(FILE *fid);
};

#endif
