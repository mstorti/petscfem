// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: autostr.h,v 1.4 2003/02/17 04:20:36 mstorti Exp $
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
  int size();
  int len();
  void resize(int m);
  void clear();
  void sprintf(char *,...);
  void vsprintf(char *,va_list ap);
  void cat_sprintf(char *,...);
  void vcat_sprintf(char *,va_list ap);
  void cat(AutoString &s);
  int getline(FILE *fid);
};

#endif
