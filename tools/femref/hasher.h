// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: hasher.h,v 1.2 2004/11/27 19:05:19 mstorti Exp $
#ifndef PETSCFEM_HASHER_H
#define PETSCFEM_HASHER_H

#include <cstdlib>

class GenHasher {
private:
  drand48_data buffer;
  unsigned short int xsubi[3];
  long int result;
public:
  GenHasher();
  void reset();
  virtual void hash(int w)=0;
  void hashv(int *w,int n);
  long int hash_val();
};

class Hasher : public GenHasher {
public:
  void hash(int w);
};

#if 0
class HasherSum {
  long int result;
public:
  void hash(int w);
};
#endif

#endif
