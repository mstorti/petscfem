// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: hasher.h,v 1.3 2004/11/27 19:26:40 mstorti Exp $
#ifndef PETSCFEM_HASHER_H
#define PETSCFEM_HASHER_H

#include <cstdlib>

class BaseHasher {
public:
  virtual void reset()=0;
  virtual void hash(int w)=0;
  virtual void hash(int *w,int n)=0;
  virtual long int val()=0;
};

class Hasher : public BaseHasher {
private:
  drand48_data buffer;
  unsigned short int xsubi[3];
  long int result;
public:
  Hasher();
  void reset();
  void hash(int w);
  void hash(int *w,int n);
  long int val();
};

class SumHasher : public BaseHasher {
private:
  drand48_data buffer;
  unsigned short int xsubi[3];
  long int result;
public:
  SumHasher();
  void reset();
  void hash(int w);
  void hash(int *w,int n);
  long int val();
};

#endif
