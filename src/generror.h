// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: generror.h,v 1.2 2003/12/03 11:52:20 mstorti Exp $
#ifndef PETSCFEM_GENERROR_H
#define PETSCFEM_GENERROR_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class GenericError : public string { 
public:
  GenericError() : string("") { }
  // GenericError(char *s) : string(s) { }
  GenericError(char *s,va_list l);
  GenericError(char *s,...);
};

#endif
