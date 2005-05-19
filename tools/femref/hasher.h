// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: hasher.h,v 1.11 2005/05/19 00:46:54 mstorti Exp $
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

class FastSumHasher : public BaseHasher {
private:
  int retval;
public:
  static int hashf(int x) {
    const unsigned int 
      c = 0x238e1f29,
      m = 0x6b8b4567;
    unsigned long int y, tmp;
#if 1
    tmp = x + c;
    y = tmp*tmp;
    y ^= m;
#else
    tmp = x ^ m;
    y = tmp*tmp;
#endif
    return y;
  }
  FastSumHasher() : retval(0) { }
  void reset() { retval=0; }
  void hash(int w) { retval += hashf(w); }
  void hash(int *w,int n) {
    for (int j=0; j<n; j++) 
      retval += hashf(w[j]); 
  }
  long int val() { return retval; }
};

class FastHasher : public BaseHasher {
private:
  const unsigned int c,m,n;
  int state;
  void hashf(int x) {
    unsigned long int tmp;
    state ^= x;
    for (int j=0; j<n; j++) {
      tmp = state + c;
      state = tmp*tmp;
      state ^= m;
    }
  }
public:
  FastHasher() 
    : state(0), 
#if 1
      c(0x238e1f29), m(0x6b8b4567),
#else
      c(54841), m(0x5c32b3a3), 
#endif
      n(1) { }
  void reset() { state = c; }
  void hash(int w) { hashf(w); }
  void hash(int *w,int n) {
    for (int j=0; j<n; j++) hashf(w[j]); 
  }
  long int val() { return state; }
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
    MD5Final(&ctx);
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

class BJHasher : public BaseHasher {
private:
  int state;
public:
  BJHasher() { reset(); }
  void reset() { state=0; }
  void hash(int w);
  void hash(int *w,int n);
  long int val();
};

#endif
