// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: linkgraph.h,v 1.3 2002/07/22 19:08:37 mstorti Exp $
#ifndef LINKGRAPH_H
#define LINKGRAPH_H

extern int MY_RANK,SIZE;

#include <src/graph.h>
#include <src/iisdgraph.h>
#include <src/dvector.h>
#include <src/graphdv.h> // for int_pair

/** Graph adjacency class that stores the list of vertices 
    adjacent to a vertex #i# as a linked list of cells pointed 
    by cursors in a dynamic vector (dvector). 
*/
class link_graph : public StoreGraph {
protected:
  /// Number of nodes in the graph
  int M;
  /// Storage array
  dvector<int_pair> da;
  /// #chunk_size# default value
  static int CHUNK_SIZE_DEF;
  /// terminator
  static int null;
  /// append to a list
  void list_insert(int header, int val);
  /// return a cursor to an available cell
  int available();
public:
  /// cuasi constructor
  void init(int MM);
  /// Constructor includes partitioner and communicator
  link_graph(int MM, const DofPartitioner *part,
	      MPI_Comm comm,int chunk_size=CHUNK_SIZE_DEF);
  /** Set new chunk size (container must be empty).
      @param new_chunk_size (input) new chunk size for the vector.
  */
  void set_chunk_size(int new_chunk_size) { da.set_chunk_size(new_chunk_size); }
  /// Destructor. Clear dynamic array. 
  ~link_graph() { da.clear(); }
  /** Add an edge to the adjacency table
    @param j (input) first index of edge
    @param k (input) second index of edge
  */
  void add(int j,int k) { list_insert(j,k); }
  /** Call back function that defines the graph. For a given
      vertex it returns the list of adjacent vertices. 
      @param v (input) the source vertex 
      @param ngbrs (output) the set of adjacent vertices 
  */
  void set_ngbrs(int v,GSet &ngbrs);
  /// Return number of edges
  int size() { return da.size(); } // not implemented yet
  /// Clear all edges
  void clear() { }
  /// Scatter among processors
  void scatter() { assert(SIZE==1); }
};

#endif
