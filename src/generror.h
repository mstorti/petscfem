// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: generror.h,v 1.1 2003/09/13 19:09:08 mstorti Exp $
#ifndef PETSCFEM_GENERROR_H
#define PETSCFEM_GENERROR_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class GenericError : public string { 
public:
  GenericError(char *s) : string(s) { }
};

#endif
