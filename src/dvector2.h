// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: dvector2.h,v 1.26 2005/01/15 23:41:32 mstorti Exp $
#ifndef PETSCFEM_DVECTOR2_H
#define PETSCFEM_DVECTOR2_H

#include <cstdarg>
#include <cassert>
#include <vector>
#define CHUNK_SIZE_INIT 10000

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>::~dvector() {
  clear();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::reff(int j,int &chunk,int &k) const {
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
dvector<T>::dvector(int cs) { 
  chunk_size = cs;
  chunks = NULL;
  size_m = nchunks = 0;
  shape_p = &*shape.begin();
  rank_m = 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::set_chunk_size(int new_chunk_size) {
  int old_size = size();
  T *tempo=NULL;
  if (old_size!=0) {
    tempo = new T[old_size];
    export_vals(tempo);
  }
  clear();
  chunk_size = new_chunk_size;
  resize(old_size);
  if (old_size!=0) {
    set(tempo);
    delete[] tempo;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T> & dvector<T>::defrag() {
  if (chunk_size<size()) set_chunk_size(size());
  assert(chunks_n()==1 || size()==0);
  return *this;
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
inline const T &dvector<T>::ref(int j) const {
  assert(j>=0 && j<size_m);
  int chunk,k;
  reff(j,chunk,k);
  return chunks[chunk][k];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::size(void) const { return size_m; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::push(const T &t) {
  int chunk,k;
  // find chunk and position
  reff(size_m,chunk,k);
  if (k==0) {
    // allocate new chunk
    assert(nchunks == chunk_vector.size());
    // resync `chunk_vector' and `nchunks'
    T* p = new T[chunk_size];
    chunk_vector.push_back(p);
    chunks = &*chunk_vector.begin();
    nchunks++;
  }
  // copy object
  ref(size_m++) = t;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::
shrink(int new_size) {
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
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>&  dvector<T>::
resize(int new_size,T &t) {
  // If new_size is greater use #push# else use #shrink#
  if (new_size > size_m) {
    while (size_m<new_size) push(t);
  } else shrink(new_size);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::
resize(int new_size) { 
  T t; 
  resize(new_size,t); 
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::
clear(void) { 
  shrink(0); 
  return *this; 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::mono(int s) {
  clear();
  set_chunk_size(s);
  resize(s);
  defrag();
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::mono(int s,T e) {
  mono(s);
  set(e);
  return *this;
}

int dvector<double>::read(FILE *fid,double &t);

int dvector<int>::read(FILE *fid,int &t);

int dvector<float>::read(FILE *fid,float &t);

int dvector<double>::printe(FILE *fid,double t);

int dvector<int>::printe(FILE *fid,int t);

int dvector<float>::printe(FILE *fid,float t);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::read(FILE *fid,T &t) {
  printf("dvector<>: error: not implemented read "
	 "function for this scalar type");
  abort();
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::printe(FILE *fid,T t) { return 1; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::read(FILE *fid) {
  // Currently for vectors only (reshape after)
  // assert(rank_m==1);
  int m=size();
  for (int j=0; j<m; j++) {
    int ierr = read(fid,ref(j));
    assert(!ierr);
  }
  if (shape.size()) 
    recompute_shape();
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::read(const char *name) {
  FILE *fid = fopen(name,"r");
  if (!fid) {
    printf("dvector<T>::read(): can't open file \"%s\"\n",name);
    abort();
  }
  read(fid);
  fclose(fid);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::cat(FILE *fid) {
  // Currently for vectors only (reshape after)
  nread = 0;
  while(1) {
    T val;
    int ierr = read(fid,val);
    if(ierr) break;
    push(val);
    nread++;
  }
  if (shape.size()) 
    recompute_shape();
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::cat(const char *name) {
  FILE *fid = fopen(name,"r");
  if (!fid) {
    printf("dvector<T>::read(): can't open file \"%s\"\n",name);
    abort();
  }
  cat(fid);
  fclose(fid);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::cat(const char *name, int &nread_a) {
  cat(name);
  nread_a = nread;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::cat(const T* in,T term) {
  const T *q = in;
  while (*q!=term) push(*q++);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::
print(FILE *fid,int rowsz) {
  int M=size();
  if (!rowsz && rank()>1) size(rank()-1);
  if (!rowsz) {
    for (int j=0; j<M; j++) {
      ierr = printe(fid,ref(j));
      if (ierr) break;
      fprintf(fid,"\n");
    }
  } else {
    int j=0;
    while (j<M) {
      ierr = printe(fid,ref(j));
      if (ierr) break;
      fprintf(fid," ");
      j++;
      if (j % rowsz == 0)
	fprintf(fid,"\n");
    }
  }
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T>& dvector<T>::
print(const char *name,int rowsz) {
  ierr = 0;
  FILE *fid = fopen(name,"w");
  if (!fid) {
    printf("dvector<T>::read(): can't open file \"%s\"\n",name);
    abort();
  }
  print(fid,rowsz);
  fclose(fid);
  assert(!ierr);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::bsearch(const T &t,int first, int last) {
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
int dvector<T>::find(const T &t,int first, int last) {
  int lastt;
  lastt = (last==-1 ? size() : last);
  int indx = bsearch(t,first,lastt);
  return indx!=lastt && ref(indx)==t;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T> &
dvector<T>::sort(int first, int last) {
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
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::remove_unique(int first, int last) {
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
dvector<T> &
dvector<T>::reshapev(int rank_a,va_list ap) {
  rank_m = rank_a;
  assert(rank_m>0);
  int new_size=1;
  shape.clear();
  for (int j=0; j<rank_m; j++) {
    int d = va_arg (ap, int);
    assert(j==0 || d>=0);
    new_size *= d;
    shape.push_back(d);
  }
  shape_p = &*shape.begin();
  assert(new_size<=0 || new_size==size());
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T> &
dvector<T>::reshape(int rank_a,...) {
  va_list ap;
  va_start(ap,rank_a);
  reshapev(rank_a,ap);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
void dvector<T>::recompute_shape() {
  int odim = 1;
  for (int k=1; k<rank_m; k++) odim *= shape[k];
  assert(size() % odim ==0);
  shape[0] = size()/odim;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::size(int dim) const {
  assert(dim<shape.size());
  return shape[dim];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::rank() const {
  return shape.size();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T> &
dvector<T>::a_resize(int rank_a,...) {
  rank_m = rank_a;
  va_list ap;
  va_start(ap,rank_a);
  assert(rank_m>0);
  int new_size=1;
  for (int j=0; j<rank_m; j++) {
    int d = va_arg (ap, int);
    new_size *= d;
    assert(d>=0);
  }
  resize(new_size);

  va_start(ap,rank_a);
  reshapev(rank_m,ap);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
const T& dvector<T>::ev(int j,va_list ap) 
  const { return ref(pos(j,ap)); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
T& dvector<T>::ev(int j,va_list ap) { 
  return ref(pos(j,ap)); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::pos(int j,va_list ap) const {
  int poss = j;
  for (int k=1; k<rank_m; k++) 
    poss = poss*shape_p[k]+va_arg(ap,int);
  return poss;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
T& dvector<T>::e(int j,...) {
  va_list ap;
  va_start(ap,j);
  return ev(j,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
const T& dvector<T>::e(int j,...) const {
  va_list ap;
  va_start(ap,j);
  return ev(j,ap);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
int dvector<T>::chunks_n() const { return nchunks; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
const T *dvector<T>::buff() const { 
  return (size() ? &ref(0) : NULL); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
T *dvector<T>::buff() { 
  return (size() ? &ref(0) : NULL); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T> &
dvector<T>::set(T t) {
  int s = size();
  for (int j=0; j<s; j++) ref(j) = t;
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
dvector<T> &
dvector<T>::set(const T* array) {
  int s = size();
  for (int j=0; j<s; j++) ref(j) = array[j];
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
template<class T>
const dvector<T> &
dvector<T>::export_vals(T* array) const {
  int s = size();
  for (int j=0; j<s; j++) array[j] = ref(j);
  return *this;
}

#endif

