// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: fastlib.h,v 1.7 2005/01/19 03:02:41 mstorti Exp $

#ifndef FASTLIB_H
#define FASTLIB_H

#define FASTLIB_CHUNK_SIZE 10
//#define TMPLARGSDEF <class T,int chunk_size=FASTLIB_CHUNK_SIZE> 

template <class T,int chunk_size=FASTLIB_CHUNK_SIZE>
class FastVector {
 public:
  inline FastVector(void) {size_=0; storage=chunk_size; 
  store=rigid_store; flexible_store=0;}; 
  inline FastVector(const int m,const T n);
  ~FastVector();
  void print(const char *s=NULL) const;
  T & operator[] (const int j) {return store[j];};
  const T & operator[] (const int j) const {return store[j];};
  inline int operator== (const FastVector & indx) const;
  inline int operator!= (const FastVector & indx) const;
  inline FastVector & operator= (const FastVector & indx);
  int size(void) const {return size_;};
  int push_back(const T j) {resize_(size_+1); store[size_++]=j; return 0;};
  void reset() {size_=0;};
  void resize(const int new_size=0) {resize_(new_size); size_=new_size;};
 private:
  inline void resize_(const int n);
  T rigid_store[chunk_size];
  T *flexible_store;
  T *store;
  int storage;
  int size_;
};

#endif
