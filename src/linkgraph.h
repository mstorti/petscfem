// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: linkgraph.h,v 1.1 2002/07/22 15:06:51 mstorti Exp $
#ifndef GRAPHDV_H
#define GRAPHDV_H

extern int MY_RANK,SIZE;

#include <src/iisdgraph.h>
#include <src/dvector.h>
#include <src/graphdv.h> // for int_pair

/** Graph adjacency class that stores the list of vertices 
    adjacent to a vertex #i# as a linked list of cells pointed 
    by cursors in a dynamic vector (dvector). 
*/
class link_graph : public StoreGraph {
protected:
  /// Storage array
  dvector<int_pair> da;
  /// #chunk_size# default value
  static int CHUNK_SIZE_DEF;
public:
  /** Constructor. May define the chunk size. 
    @author M. Storti
    @param chunk_size (input) size of chunks used in the
    internal dynamic vector
  */
  graphdv(int chunk_size=CHUNK_SIZE_DEF);
  /** Set new chunk size (container must be empty).
      @param new_chunk_size (input) new chunk size for the vector.
  */
  void set_chunk_size(int new_chunk_size) { da.set_chunk_size(new_chunk_size); }
  /// Destructor. Clear dynamic array. 
  ~graphdv() { da.clear(); }
  /** Add an edge to the adjacency table
    @param j (input) first index of edge
    @param k (input) second index of edge
  */
  void add(int j,int k);
  /** Prints adjacency table (first resync)
      @param s (input) A string to be printed along with the adjacency graph
  */
  void print(const char *s=NULL);
  /** Call back function that defines the graph. For a given
      vertex it returns the list of adjacent vertices. 
      @author M. Storti
      @param v (input) the source vertex 
      @param ngbrs (output) the set of adjacent vertices 
  */
  void set_ngbrs(int v,vector<int> &ngbrs);
  /** Call back function that defines the graph. For a given
      vertex it returns the list of adjacent vertices. 
      @author M. Storti
      @param v (input) the source vertex 
      @param ngbrs (output) the set of adjacent vertices 
  */
  void set_ngbrs(int v,GSet &ngbrs);
  /// Return number of edges
  int size() { resync(); return da.size(); }
  /// Clear all edges
  void clear() { da.clear(); modif=0; ordered=0; max = MAX_INIT; }
  /// Constructor includes partitioner and communicator
  link_graph(int M, const DofPartitioner *part,
	      MPI_Comm comm,int chunk_size=CHUNK_SIZE_DEF) : graphdv(chunk_size) {}
  /// Scatter among processors
  void scatter() { assert(SIZE==1); }
};

#endif
