// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: dvector.h,v 1.1 2002/07/21 13:02:14 mstorti Exp $
#ifndef DVECTOR_H
#define DVECTOR_H

template<class T>
class dvector {
private:
  /// size of each chunk
  int chunk_size;
  /// total number of elements 
  int size_m;
  /// Number of chunks position of first free element
  int nchunks;
  /// pointers to chunks
  vector<T *> chunk_vector;
  /// pointer to STL storage array in chunk_vector
  T **chunks;
  /** Find in what chunk, in what position 'k' is 'j'
      @param j (input) the position in the vector
      @param chunk (output) number of chunk
      @param k (output) the position in chunk
  */ 
  void reff(int j,int &chunk,int &k) {
    chunk = j/chunk_size;
    k = j % chunk_size;
  }
  /** Exchange objescts `u' and `v'
      @param u (input/output) first element to be exchanged
      @param v (input/output) second element to be exchanged
  */ 
  void X(T &u,T &v) {
    T t = u;
    u = v;
    v = t;
  }
  /** The origial pascal code for the heap-sort was written in pascal
      with 1-based indices.. The macro #V# shifts to map index 1 to #first#
      @param j (input) base 1 index 
      @return 1 to #first# mapped index 
  */ 
#define V(j) ref(first+(j)-1)
  /**  Auxiliary function to the ``sort'' (heap-sort) routine. 
       Note: The heap is maximal. 
      @param first (input) first position of the heap
      @param p (input) position to re-heap
      @param u (input) last position in heap
  */ 
  void push_heap (int first, int p, int u) {
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
public:
  /** Constructor 
      @param cs (input) chunk size for this vector
      @return a reference to the matrix.
  */  
  dvector(int cs=100) { 
    chunk_size = cs;
    chunks = NULL;
    size_m = nchunks = 0;
  }
  /** reference to the object in some position
      @param j (input) position in the vector
      @return a reference to the object inside the vector 
  */ 
  T &ref(int j) {
    assert(j>=0 && j<size_m);
    int chunk,k;
    reff(j,chunk,k);
    return chunks[chunk][k];
  }
  /** number of elements in the vector
      @return size of vector
  */ 
  int size(void) { return size_m; }
  /** Insert an element at the last position
      @param t (input) element to be appended
  */ 
  void push(const T &t) {
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
  /** Resizes but only for a new dimension 
      lower than the original. This can be used even if there is not a
      default constructor for the base class. 
      @param new_size (input) new size for the vector
  */
  void shrink(int new_size) {
    assert(new_size<=size_m);
    int chunk,k;
    // demangle new size
    reff(new_size,chunk,k);
    // free all chunks above the last needed
    while (nchunks > chunk+1) {
      delete[] chunk_vector[nchunks-1];
      chunk_vector.pop_back();
      nchunks--;
    }
    size_m = new_size;
  }
  /** Resizes the vector.
      @param new_size (input) new size for the vector
  */
  void resize(int new_size,T t) {
    // If new_size is greater use #push# else use #shrink#
    if (new_size > size_m) {
      while (size_m<new_size) push(t);
    } else shrink(new_size);
  }
  /** Resizes the vector with default constructor.
      @param new_size (input) new size for the vector
  */ 
  void resize(int new_size) { resize(new_size,T()); }
  /// Resizes to null size.
  void clear(void) { shrink(0); }
  /**  Sort the elements in range #first#, #last-1#. 
       Uses ``heap-sort'' algorithm. 
       @param first (input) first position in range to be sorted
       @param last (input) past to the end position in range to be sorted
  */  
  int (const T &t,int first=0, int last=-1) {
  /** Binary search algorithm. Finds #t# in the range #first# to
      #last-1#.  
      @param t (input) element to find
      @param first (input) first position in range to search
      @param last (input) past to the end  position in range to search
      @return position of first element greater of equal than the
      searched element.
  */ 
  int bsearch(const T &t,int first=0, int last=-1) {
    if (last==-1) last=size();
    if (ref(first)>=t) return first;
    if (last<=first) return last;
    int p=first, q=last, r;
    while (1) {
      if (q==p+1) return q;
      r = (p+q)/2;
      if (ref(r)<t) p=r;
      else q=r;
    }
    assert(0);
  }
  /** Finds whether a given element is or not in the vector.
      @param t (input) element to find
      @param first (input) first position in range to search
      @param last (input) past to the end  position in range to search
      @return true/false if the element is or is not in the vector. 
  */ 
  int find(const T &t,int first=0, int last=-1) {
    return find(t,first,last)!=last;
  }
  /** 
      @param (input)
      @return a reference to the matrix.
  */ 
  void sort(int first=0, int last=-1) {
    if (last==-1) last=size();
    int i, j, n = last-first;
    T t;
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
};

#endif
