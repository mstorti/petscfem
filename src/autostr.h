// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: autostr.h,v 1.1 2003/02/16 03:15:26 mstorti Exp $
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
  void cat(AutoString &s);
};

#endif
