// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: hasher.h,v 1.7 2004/11/28 23:33:32 mstorti Exp $
#ifndef PETSCFEM_HASHER_H
#define PETSCFEM_HASHER_H

#include <cstdlib>
extern "C" {
#include "./md5.h"
}

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

class PowSumHasher : public BaseHasher {
private:
  int retval;
public:
  static int hashf(int x) {
    const unsigned int 
      c = 0x238e1f29,
      m = 0x6b8b4567,
      n = 0;
    unsigned long int y, tmp;
    y = x;
    for (int j=0; j<n; j++) {
      tmp = y + c;
      y = tmp*tmp;
      y ^= m;
    }
    return y;
  }
  PowSumHasher() : retval(0) { }
  void reset() { retval=0; }
  void hash(int w) { retval += hashf(w); }
  void hash(int *w,int n) {
    for (int j=0; j<n; j++) 
      retval += hashf(w[j]); 
  }
  long int val() { return retval; }
};

class BasicSumHasher : public BaseHasher {
private:
  int retval;
public:
  BasicSumHasher() : retval(0) { }
  void reset() { retval=0; }
  void hash(int w) { retval += w; }
  void hash(int *w,int n) {
    for (int j=0; j<n; j++) 
      retval += w[j]; 
  }
  long int val() { return retval; }
};

class MD5Hasher : public BaseHasher {
private:
  MD5_CTX ctx;
public:
  MD5Hasher() { reset(); }
  void reset() { MD5Init(&ctx); }
  void hash(int w) { 
    MD5Update(&ctx,(unsigned char *)&w,sizeof(int)); 
  }
  void hash(int *w,int n) { 
    MD5Update(&ctx,(unsigned char *)w,n*sizeof(int)); 
  }
  long int val() {
    int retval;
    memcpy(&retval,ctx.digest,sizeof(int));
    return retval;
  }
};

class MD5SumHasher : public BaseHasher {
private:
  MD5_CTX ctx;
  int retval;
public:
  MD5SumHasher() { reset(); }
  void reset() { retval=0; }
  void hash(int w) { 
    MD5Init(&ctx);
    MD5Update(&ctx,(unsigned char *)&w,sizeof(int)); 
    MD5Final(&ctx);
    int h;
    memcpy(&h,ctx.digest,sizeof(int));
    retval += h;
  }
  void hash(int *w,int n) { 
    for (int j=0; j<n; j++)
      hash(w[j]);
  }
  long int val() {
    return retval;
  }
};

#endif
