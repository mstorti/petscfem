/*__INSERT_LICENSE__*/
// $Id: tryme3.cpp,v 1.6 2002/07/20 17:24:06 mstorti Exp $

#include <cmath>
#include <map>
#include <set>
#include <vector>

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

template < class T > class testalloc;

template <> class testalloc< void > {
public:
  typedef void* pointer;
  typedef const void* const_pointer;
  typedef void value_type;
  template < class U > struct rebind { typedef testalloc<U> other; };
};

//#define TESTALLOC_DEBUG

template < class T > class testalloc {
private:
  
public:
  typedef T value_type;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef T& reference;
  typedef const T& const_reference;
  pointer address(reference r) const {
#ifdef TESTALLOC_DEBUG
    cout << "testalloc<T>::address called with r = " 
	 << r << "." << endl;
#endif
    return &r;
  }
  const_pointer address( const_reference r ) const {
#ifdef TESTALLOC_DEBUG
    cout << "testalloc<T>:address const called with r = " << r 
	 << "." << endl;
#endif
    return &r;
  }
  testalloc( void ) throw() {
#ifdef TESTALLOC_DEBUG
    cout << "testalloc(void) ctor being called." << endl;
#endif
  }
  template < class U > testalloc( const testalloc<U>& ) throw() {
#ifdef TESTALLOC_DEBUG
    cout << "testalloc() copy-ctor being called." << endl;
#endif
  }
  ~testalloc( void ) throw() {
#ifdef TESTALLOC_DEBUG
    cout << "~testalloc dtor being called." << endl;
#endif
  }
  pointer allocate( size_type n, testalloc<void>::const_pointer hint = 0 ) {
#ifdef TESTALLOC_DEBUG
    cout << "testalloc<T>::allocate( size_type = " << n << ", hint ="
	 << hint << " )." << endl;
#endif
    if( n>1 ) {
      return new T[n];
    } else {
      return new T();
    }
  }
  void deallocate( pointer p, size_type n ) {
#ifdef TESTALLOC_DEBUG
    cout << "testalloc<T>::deallocate( pointer = " << p
	 << ", size_type = " << n << " ) being called." << endl;
#endif
#if 0
    if( p ) {
      if( n>1 ) {
	delete[](p);
      } else {
	delete(p);
      }
    }
#endif
  }
  void construct( pointer p, const T& val ) {
#ifdef TESTALLOC_DEBUG
    cout << "testalloc<T>::construct( pointer = " << p
	 << ", const T& = " << val << " ) being called." << endl;
#endif
    new(p) T(val);
  }
  void destroy( pointer p ) {
#ifdef TESTALLOC_DEBUG
    cout << "testalloc<T>::destroy( pointer = " << p 
	 << " ) being called." << endl;
#endif
    p->~T();
  }
  size_type max_size( void ) const throw() { 
    return ((size_t)-1)/sizeof( T ); 
  }
  template < class U > struct rebind { typedef testalloc<U> other; };
};

//#define ALLOC malloc_alloc
//#define ALLOC alloc
#define ALLOC testalloc<int>
typedef set< int, less<int>, ALLOC > SET;

class MAP : public map< int, SET, less<int>, ALLOC > {
public:
  void add(int i,int j) { (*this)[i].insert(j); }
};

int main(int argc, char **argv) {
  int M=10000;
  MAP g;
  for (int j=0; j<M; j++) { g.add(irand(1,M),irand(1,M)); }
  double *a = new double[3*M];
  for (int j=0; j<3*M; j++) a[j] = double(j);
  g.clear();
  delete[] a;
}
