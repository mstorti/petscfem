// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: fastlib2.h,v 1.3 2003/07/03 04:32:11 mstorti Exp $
#ifndef FASTLIB2_H
#define FASTLIB2_H

using namespace std;

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
inline int FastVector <T,chunk_size>::operator!= (const FastVector & indx) const {
  return !(*this == indx);
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

#endif
