// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: hasher.h,v 1.1 2004/11/24 23:22:35 mstorti Exp $
#ifndef PETSCFEM_HASHER_H
#define PETSCFEM_HASHER_H

#include <cstdlib>

class Hasher {
private:
  drand48_data buffer;
  unsigned short int xsubi[3];
  long int result;
public:
  Hasher();
  void reset();
  void hash(int w);
  void hash(int *w,int n);
  long int hash_val();
};

#endif
