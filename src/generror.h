// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: generror.h,v 1.4 2003/12/06 15:11:09 mstorti Exp $
#ifndef PETSCFEM_GENERROR_H
#define PETSCFEM_GENERROR_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class GenericError : public string { 
public:
  GenericError() : string("") { }
  GenericError(string s) : string(s) { }
  // GenericError(char *s) : string(s) { }
  GenericError(char *s,va_list l);
  GenericError(char *s,...);
};

#endif
