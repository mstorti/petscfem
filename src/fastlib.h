// -*- mode: c++ -*-
//__INSERT_LICENSE__
//$Id: fastlib.h,v 1.3 2001/04/01 01:35:06 mstorti Exp $

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
  inline ~FastVector();
  void print(const char *s=NULL) const;
  T & operator[] (const int j) {return store[j];};
  const T & operator[] (const int j) const {return store[j];};
  inline int operator== (const FastVector & indx) const;
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template <class T,int chunk_size>
FastVector <T,chunk_size> & FastVector <T,chunk_size> ::operator= (const FastVector & indx) {
  store = rigid_store;
  if (flexible_store) {
    delete[] flexible_store;
  }
  resize_(indx.size_);
  size_ = indx.size_;
  for (int j=0; j<size_; j++) {
    store[j] = indx.store[j];
  }
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template <class T,int chunk_size>
inline void FastVector <T,chunk_size>::resize_(const int n) {
  if (n <= storage) {
    return;
  } else {
    storage += chunk_size;
    if (storage<n) storage = n;
    T *new_flexible_store = new T[storage];
    for (int j=0; j<size_; j++) 
      new_flexible_store[j] = store[j];
    if (flexible_store) {
      delete[] flexible_store;
    }
    store = flexible_store = new_flexible_store;
  }
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template <class T,int chunk_size>
inline FastVector <T,chunk_size>::~FastVector() {
  if(flexible_store) {
    delete[] flexible_store;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template <class T,int chunk_size>
inline int FastVector <T,chunk_size>::operator== (const FastVector & indx) const {
  if (size_!=indx.size_) return 0;
  for (int j=0; j<size_; j++) {
    if (store[j] != indx.store[j]) return 0;
  }
  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template <class T,int chunk_size>
inline FastVector <T,chunk_size>::FastVector(const int m,const T n) {
  // flexible_store=NULL;
  flexible_store=0;
  store=rigid_store;
  storage=chunk_size;
  size_=0;
  resize_(m);
  size_=m;
  for (int j=0; j<m; j++) {
    store[j]=n;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template <class T,int chunk_size>
void FastVector <T,chunk_size>::print(const char *s) const {
  if (s!=NULL) cout << s << "  ";
  int ndims = size();
  for (int jd=0; jd<ndims; jd++) {
    cout << (*this)[jd] << "  ";
  }
  cout << endl;
}

#if 0
template <class Key,class Value>
class Pair {
 public:
  Key first;
  Value second;
};

template <class Key,class Value,int chunk_size=FASTLIB_CHUNK_SIZE>
class FastMap : public FastVector<Pair<Key,Value>> {
 public:
  inline FastVector(void) {size_=0; storage=chunk_size; 
  store=rigid_store; flexible_store=0;}; 
  inline FastVector(const int m,const T n);
  inline ~FastVector();
  void print(const char *s=NULL) const;
  T & operator[] (const int j) {return store[j];};
  const T & operator[] (const int j) const {return store[j];};
  inline int operator== (const FastVector & indx) const;
  inline FastVector & operator= (const FastVector & indx);
  int size(void) const {return size_;};
  int push_back(const T j) {resize_(size_+1); store[size_++]=j;};
  void reset() {size_=0;};
 private:
  inline void resize_(const int n);
  T rigid_store[chunk_size];
  T *flexible_store;
  T *store;
  int storage;
  int size_;
};
#endif

#endif
