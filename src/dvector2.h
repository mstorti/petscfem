// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: dvector2.h,v 1.3 2003/02/25 20:34:22 mstorti Exp $
#ifndef PETSCFEM_DVECTOR2_H
#define PETSCFEM_DVECTOR2_H

#include <cstdarg>
#include <vector>
#define CHUNK_SIZE_INIT 10000

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::reff(int j,int &chunk,int &k) {
  chunk = j/chunk_size;
  k = j % chunk_size;
  if (k<0) {
    chunk -=1;
    k += chunk_size;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::X(T &u,T &v) {
  T t = u;
  u = v;
  v = t;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** The original pascal code for the heap-sort was written in pascal
    with 1-based indices.. The macro #V# shifts to map index 1 to #first#
    @param j (input) base 1 index 
    @return 1 to #first# mapped index 
*/ 
#define V(j) ref(first+(j)-1)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::push_heap (int first, int p, int u) {
  // i1,i2,q,r are all positions in the heap (base 1)
  // the macro `V' shifts to map 1 -> first
  int i1, i2, q, r;
  r = p ; // {indica posicion actual de V [primero] }
  q = u / 2;
  while (r <= q ) {
    // i1,i2 ares left/right sons of r
    i1 = 2 * r ;
    i2 = 2 * r + 1 ;
    if (u == i1) {
      // r has only one son
      if (V(r) < V(i1)) X( V(r), V(i1));
      break;
    } else if ( V(r) < V(i1) && V(i1) >= V(i2)) {
      // r descends left
      X( V(r), V(i1)) ;
      r = i1 ;
    } else if ( V(r) < V(i2) && V(i2) >  V(i1)) {
      // r descends right
      X( V(r), V(i2)) ;
      r = i2 ;
    } else break; // Stop descending
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>::dvector(int cs = CHUNK_SIZE_INIT) { 
  chunk_size = cs;
  chunks = NULL;
  size_m = nchunks = 0;
  shape_p = shape.begin();
  rank = 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::set_chunk_size(int new_chunk_size) {
  int old_size = size();
  T *tempo;
  if (old_size!=0) {
    tempo = new T[old_size];
    for (int j=0; j<old_size; j++) tempo[j] = ref(j);
  }
  clear();
  chunk_size = new_chunk_size;
  resize(old_size);
  if (old_size!=0) {
    for (int j=0; j<old_size; j++) ref(j) = tempo[j];
  }
  delete[] tempo;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::defrag() {
  if (chunk_size<size()) set_chunk_size(size());
  assert(chunks_n()==1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
inline T &dvector<T>::ref(int j) {
  assert(j>=0 && j<size_m);
  int chunk,k;
  reff(j,chunk,k);
  return chunks[chunk][k];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::size(void) { return size_m; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::push(const T &t) {
  int chunk,k;
  // find chunk and position
  reff(size_m,chunk,k);
  if (k==0) {
    // allocate new chunk
    assert(nchunks == chunk_vector.size());
    // resync `chunk_vector' and `nchunks'
    chunk_vector.push_back(new T[chunk_size]);
    chunks = chunk_vector.begin();
    nchunks++;
  }
  // copy object
  ref(size_m++) = t;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::shrink(int new_size) {
  assert(new_size<=size_m);
  int chunk,k;
  // demangle new size
  reff(new_size-1,chunk,k);
  // free all chunks above the last needed
  while (nchunks > chunk+1) {
    delete[] chunk_vector[nchunks-1];
    chunk_vector.pop_back();
    nchunks--;
  }
  size_m = new_size;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::resize(int new_size,T &t) {
  // If new_size is greater use #push# else use #shrink#
  if (new_size > size_m) {
    while (size_m<new_size) push(t);
  } else shrink(new_size);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::resize(int new_size) { T t; resize(new_size,t); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::clear(void) { shrink(0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::bsearch(const T &t,int first=0, int last=-1) {
  int lastt=last;
  if (lastt==-1) lastt=size();
  if (lastt<=first) return lastt; // vector is void
  if (ref(first)>=t) return first;
  int p=first, q=lastt, r;
  while (1) {
    if (q==p+1) return q;
    r = (p+q)/2;
    if (ref(r)<t) p=r;
    else q=r;
  }
  assert(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::find(const T &t,int first=0, int last=-1) {
  int lastt;
  lastt = (last==-1 ? size() : last);
  int indx = bsearch(t,first,lastt);
  return indx!=lastt && ref(indx)==t;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::sort(int first=0, int last=-1) {
  if (last==-1) last=size();
  int i, j, n = last-first;
  j = n / 2 ;
  // make heap (from bottom to root)
  for (i = j; i>=1; i--) push_heap(first,i,n) ;
  // extracts root of heap (the maximum) puts at the back
  // and reheap
  for (i = n; i>=2;  i--) {
    X(V(1),V(i));
    push_heap(first,1,i-1);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::remove_unique(int first=0, int last=-1) {
  if (last==-1) last=size();
  sort(first,last);
  int p=first, e = last, q;
  if (p==e) return first;	// vector is void
  q=p;
  // move pointer `q' and when finding a differnt element copy to a
  // new position at `p'
  while (++q!=e) 
    if (ref(q)!=ref(p)) ref(++p) = ref(q);
  return p+1;
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::reshape(int rank_a,va_list ap) {
  rank = rank_a;
  assert(rank>0);
  int new_size=1;
  shape.clear();
  for (int j=0; j<rank; j++) {
    int d = va_arg (ap, int);
    assert(j==0 || d>0);
    new_size *= d;
    shape.push_back(d);
  }
  shape_p = shape.begin();
  assert(new_size<0 || new_size==size());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::reshape(int rank_a,...) {
  va_list ap;
  va_start(ap,rank_a);
  reshape(rank_a,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::a_resize(int rank_a,...) {
  rank = rank_a;
  va_list ap;
  va_start(ap,rank_a);
  assert(rank>0);
  int new_size=1;
  for (int j=0; j<rank; j++) {
    int d = va_arg (ap, int);
    new_size *= d;
    assert(d>0);
  }
  resize(new_size);

  va_start(ap,rank_a);
  reshape(rank,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
const T& dvector<T>::e(int j,va_list ap) const { return ref(0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
T& dvector<T>::e(int j,va_list ap) {
  int pos = j;
  for (int j=1; j<rank; j++) pos = pos*shape_p[j]+va_arg(ap,int);
  return ref(pos);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
T& dvector<T>::e(int j,...) {
  va_list ap;
  va_start(ap,j);
  return e(j,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::chunks_n() const { return nchunks; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
const T *dvector<T>::buff() const { return &ref(0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
T *dvector<T>::buff() { return &ref(0); }

#endif
