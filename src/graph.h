// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: graph.h,v 1.17.42.1 2005/02/18 09:02:47 mstorti Exp $
#ifndef GRAPH_H
#define GRAPH_H

#include <queue>
#include <set>
#include <memory>

using namespace std;

// Allocators, from http://www.sgi.com

// alloc: The default allocator. It is thread-safe, and usually has
// the best performance characteristics.

// pthread_alloc: A thread-safe allocator that uses a different memory
// pool for each thread; you can only use pthread_alloc if your
// operating system provides pthreads. Pthread_alloc is usually faster
// than alloc, especially on multiprocessor systems. It can, however,
// cause resource fragmentation: memory deallocated in one thread is
// not available for use by other threads.

// single_client_alloc: A fast but thread-unsafe allocator. In
// programs that only have one thread, this allocator might be faster
// than alloc.

//  malloc_alloc: An allocator that simply uses the standard library
//  function malloc. It is thread-safe but slow; the main reason why
//  you might sometimes want to use it is to get more useful
//  information from bounds-checking or leak-detection tools while you
//  are debugging.

#if __GNUC__ > 2
#define STL_ALLOCATOR __single_client_alloc 
//#define STL_ALLOCATOR __alloc
//#define STL_ALLOCATOR __new_alloc
//#define STL_ALLOCATOR __malloc_alloc
#else
#define STL_ALLOCATOR alloc
#endif

/// A set of neighbors. 
//typedef set<int,less<int>,STL_ALLOCATOR> GSet;
typedef set<int> GSet;

class Graph {
  /*** ngbrs:= Auxiliary queue that stores those nodes that are waiting
       to be marked. 
  */
  queue<int> ngbrs; 
  // adjncy:= xadj:= graph described in CSR format (as defined in
  // Metis documentation)
  int *el2vrtx;
  // vpart:= partitioning of vertices
  int *vpart;
  /// The number of vertices in the fine `graph'
  int nvrtx_f;
  /// The number of vertices in the coarse `graph'
  int nvrtx;
  /// max number of vertices in the allowed in the coarse `graph'
  int max_partgraph_vertices;
  /// The partitioning in the fine mesh (if needed)
  int *vpartf;
 public:
  /// Vertex weights 
  double weight_scale;
  /// Constructor
  Graph(int N=0) { init(N); };
  /// Set number of vertices
  virtual void init(int N);
  /** Virtual destructor (Just to prevent instantiations of this
      class. See Thinking in C++, by Bruce Eckel).
  */
  virtual ~Graph()=0;
  /// partitioning routine (calls Metis)
  void part(int max_partgraph_vertices,
	    int npart,float *tpwgts=NULL);
  /// return an array with the domain indices of all fine vertices 
  const int *vrtx_part();
  /// return the partition index for a given fine graph vertex 
  int vrtx_part(int);
  /** callback user function to return the neighbors for a 
      fine vertex `vrtx_f'. #ngbrs_v# is cleared before calling
      #set_ngbrs()#. 
  */
  virtual void set_ngbrs(int vrtx_f,GSet &ngbrs_v)=0;
  /// Callback user function for the user to set the weight of a given fine vertex. 
  virtual double weight(int vrtx_f);
  /// Clean all memory related 
  virtual void clear();
  /// Print the graph
  void print();
};

#endif
