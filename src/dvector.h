// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: dvector.h,v 1.29 2005/06/06 21:49:09 mstorti Exp $
#ifndef DVECTOR_H
#define DVECTOR_H
#include <cstdarg>
#include <cstdio>
#include <vector>
#define CHUNK_SIZE_INIT 10000

using namespace std;

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
  /// Defines the shape for accessing through #[]#
  vector<int> shape;
  /// Pointer to the shape vector
  int *shape_p;
  /// rank of the tensor (number of indices)
  int rank_m;
  /// pointer to STL storage array in chunk_vector
  T **chunks;
  /** Find in what chunk, in what position 'k' is 'j'
      @param j (input) the position in the vector
      @param chunk (output) number of chunk
      @param k (output) the position in chunk */ 
  inline void reff(int j,int &chunk,int &k) const;

  /** Exchange objects `u' and `v'
      @param u (input/output) first element to be exchanged
      @param v (input/output) second element to be exchanged */ 
  void X(T &u,T &v);

  /**  Auxiliary function to the ``sort'' (heap-sort) routine. 
       Note: The heap is maximal. 
      @param first (input) first position of the heap
      @param p (input) position to re-heap
      @param u (input) last position in heap */ 
  void push_heap (int first, int p, int u);

  /// Reads an element from a file
  int read(FILE *fid,T &t);

  /// Prints an element to a stream
  int printe(FILE *fid,T t);

  /// Prints an element to a stream
  T new_t_elem();

  /// Number of elements read with cat, read and friends
  int nread;
  /// Internal error code
  int ierr;

  /** Same as #ev(i,ap)# but used interally. Both versions,
      const and non-const refer internally to this. 
      @param j (input) first index
      @param ap (input) remaining indices
      @return reference to the element */ 
  int pos(int j,va_list ap) const;

  /** Same as #pos(j,ap)# but indices are passed 
      via a #vector<int>#. 
      @param indx (input) the indices of position. 
      @return reference to the element */ 
  int pos(vector<int> &indx) const;

public:
  /** Constructor 
      @param cs (input) chunk size for this vector
      @return a reference to the matrix. */  
  dvector(int cs = CHUNK_SIZE_INIT);

  /** Set new chunk size.
      @param new_chunk_size (input) new chunk size for the vector. */
  void set_chunk_size(int new_chunk_size);

  /** Set new chunk size in order to have onely one chunk. */
  dvector<T> & defrag();

  /** Destructor. Frees memory. */
  ~dvector();

  /** reference to the object in some position
      @param j (input) position in the vector
      @return a reference to the object inside the vector */ 
  T &ref(int j);

  /** reference to the object in some position (const version)
      @param j (input) position in the vector
      @return a reference to the object inside the vector */ 
  const T &ref(int j) const;

  /** number of elements in the vector
      @return size of vector */ 
  int size(void) const;

  /** Gives number of chunks, in which internally is fragmented 
      the data. May be used in order to check whether the
      buffer may be used as a contiguous area of memory. 
      @return number of chunks. Should be 1 if the buffer is not fragmented */ 
  int chunks_n() const;

  /** Gives access to internal buffer. 
      @return pointer to internal data area */ 
  const T *buff() const;

  /** Gives access to internal buffer. 
      @return pointer to internal data area */ 
  T *buff();

  /** Insert an element at the last position
      @param t (input) element to be appended */ 
  dvector<T>& push(const T &t);

  /** Resizes but only for a new dimension 
      lower than the original. This can be used even if there is not a
      default constructor for the base class. 
      @param new_size (input) new size for the vector */
  dvector<T>& shrink(int new_size);

  /** Resizes the vector.
      @param new_size (input) new size for the vector */
  dvector<T>&  resize(int new_size,T &t);

  /** Resizes the vector with default constructor.
      @param new_size (input) new size for the vector */ 
  dvector<T>& resize(int new_size);

  /// Resizes to null size.
  dvector<T>& clear(void);

  /// Resizes monolithically
  dvector<T>& mono(int size);

  /// Resizes monolithically and initializes
  dvector<T>& mono(int size,T e);

  /** Read from a stream, like Octave's `aload'.
      @param fid (input) stream to read in
      @return reference to self  */ 
  dvector<T>& read(FILE *fid);

  /** Read from a file, given its name, like Octave's `aload'.
      @param name (input) name of file to read
      @return reference to self  */ 
  dvector<T>& read(const char *name);

  /** Read from a binary file, given its name.
      @param name (input) name of file to read
      @return reference to self  */ 
  dvector<T>& readb(const char *name);

  /** Write to a binary file, given its name.
      @param name (input) name of file to write
      @return reference to self  */ 
  dvector<T>& writeb(const char *name);

  /** Read from a stream, as much as possible,
      appending to the vector */
  dvector<T>& cat(FILE *fid);

  /** Read from a file, as much as possible, 
      appending to the vector */
  dvector<T>& cat(const char *name);

  /** Read from a file, given its name, like Octave's `aload'.
      @param name (input) name of file to read
      @param nread (output) number of elments read
      @return reference to self  */ 
  dvector<T>& cat(const char *name,int &nread);

  /** Read from a standard C vector, until `term' is found. 
      @param in (input) the C vector of objects from 
      where to fill the container. 
      @param term (output) the element that indicates the end. 
      @return reference to self  */ 
  dvector<T>& cat(const T* in,T term);

  /** Clone from #w#, with the same size and shape, 
      making #*this=w#. 
      @param w (input) the vector to clone. 
      @return reference to self  */ 
  dvector<T>& clone(const dvector<T> &w);

  /** Perform the AXPY operation #*this += alpha*y#
      @param alpha (input) the scalar parameter
      @param y (input) the vector
      @return reference to self */ 
  dvector<T>& axpy(double alpha,const dvector<T> &y);

  /** Scale #*this# by #alpha#. 
      @param alpha (input) the scalar to scale. 
      @return reference to self  */ 
  dvector<T>& scale(double alpha);

  /// Write to a stream
  dvector<T>& print(FILE *fid=stdout,int rowsz=0);

  /// Write to a file by name
  dvector<T>& print(const char *name,int rowsz=0);

  /** Binary search algorithm. Finds #t# in the range #first# to
      #last-1#.  
      @param t (input) element to find
      @param first (input) first position in range to search
      @param last (input) past to the end  position in range to search
      @return position of first element greater or equal than the
      searched element. (Assume a infty at end). */ 
  int bsearch(const T &t,int first=0, int last=-1);

  /** Finds whether a given element is or not in the vector.
      @param t (input) element to find
      @param first (input) first position in range to search
      @param last (input) past to the end  position in range to search
      @return true/false if the element is or is not in the vector. */ 
  int find(const T &t,int first=0, int last=-1);

  /**  Sort the elements in range #first#, #last-1#. 
       Uses ``heap-sort'' algorithm. 
       @param first (input) first position in range to be sorted
       @param last (input) past to the end position in range to be sorted */  
  dvector<T> &sort(int first=0, int last=-1);

  /** Sorts elements in a vector range and eliminates repeated values. 
      @param first (input) first position in range to be sorted
      @param last (input) past to the end position in range to be sorted
      @return past the end position to the ordered region */ 
  int remove_unique(int first=0, int last=-1);

  /** Reshapes the tensor with #rank# dimensions, and the rest of arguments
      are the size of each dimension.
      @param rank (input) number of dimensions
      @param j,k,l (input) the sizes of each argument */
  dvector<T> & reshape(int rank,...);

  /** Reshapes the tensor with #rank# dimensions, and the rest of arguments
      are the size of each dimension.
      @param shape (input) an integer vector with the shape. */
  dvector<T> & reshape(const vector<int> &shape);

  /** Export the shape of the dvector to vector #shape#. 
      @param shape (output) External container where the shape is exported. */ 
  void get_shape(vector<int> &shape) const;

  /** Same as #reshape(i,j,k,l)# but the variable argument list 
      is pased as a variadic macro. 
      @param rank (input) number of dimensions. 
      @param ap (input) the list of dimensions, as a variadic macro */ 
  dvector<T> &reshapev(int rank_a,va_list ap);

  /** Returns the dimension of the #dim# index. 
      @param dim (input) the index. 
      @return The dimension for index #dim# */ 
  int size(int dim) const;

  /** Returns the rank of the array. 
      @return the rank of the array. */ 
  int rank() const;

  /** Recompute the shape vector. Verify that
      #size()# is compatible with the given shape. */ 
  void recompute_shape();

  /** Array resize the dynamic vector. Resizes and reshapes the vector
      by defining its rank and list of dimensions. 
      @param rank (input) number of dimensions 
      @param i,j,k (input) list of dimensions */ 
  dvector<T> &a_resize(int rank,...);

  /** Gets a particular element of the array
      @param i,j,k,l (input) indices 
      @return a reference to the internal element */ 
  T& e(int j,...);

  /** Gets a particular element of the array
      @param indx (input) vector of indices 
      @return a reference to the internal element */ 
  T& e(vector<int> &indx);

  /** Gets a particular element of the array. Const version. 
      @param indx (input) vector of indices 
      @return a reference to the internal element */ 
  const T& e(vector<int> &indx) const;

  /** Same as #e(i,j,k,...)# but the list of arguments is passed as
      a variadic macro. 
      @param j (input) first index
      @param ap (input) remaiing indices
      @return reference to the element */ 
  T& ev(int j,va_list ap);

  /** Same as #ev(i,j,k,...)# but returns a const reference. 
      @param i,j,k,l (input) indices 
      @return const reference to the internal element */ 
  const T& e(int j,...) const;

  /** Same as #ev(int,va_list)# but returns a const reference. 
      @param j (input) first index
      @param ap (input) remaiing indices
      @return const reference to the element */ 
  const T& ev(int j,va_list ap) const;

  /** Sets all values from an array.
      @param array (input) the values to be set. */ 
  dvector<T> &set(const T* array);

  /** Sets all elements to a constant value. 
      @param t (input) the value to be set. */ 
  dvector<T> & set(T t);

  /** Exports all values to an external array.
      @param array (input) the buffer where to put the values. */ 
  const dvector<T> & export_vals (T* array) const;

};

#endif
