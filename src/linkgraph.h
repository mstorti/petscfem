// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: linkgraph.h,v 1.5 2002/07/23 03:42:38 mstorti Exp $
#ifndef LINKGRAPH_H
#define LINKGRAPH_H

extern int MY_RANK,SIZE;

#include <src/distcont.h>
#include <src/graph.h>
#include <src/iisdgraph.h>
#include <src/dvector.h>
#include <src/graphdv.h> // for int_pair

/** Graph adjacency class that stores the list of vertices 
    adjacent to a vertex #i# as a linked list of cells pointed 
    by cursors in a dynamic vector (dvector). 
*/
class link_graph {
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
  class row : public set<int> {};
  // class iterator : public int {};
  typedef int iterator;
  typedef int const_iterator;
  int begin() { return 0; }
  iterator end() { return M; }
  void erase(int first,int last) {}
  /// cuasi constructor
  void init(int MM);
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
  void clear() { M=0; da.clear(); }
};

typedef link_graph::row link_graph_row; 
typedef Partitioner<link_graph_row> link_graph_row_part;

typedef  
DistCont<link_graph,link_graph_row,link_graph_row_part>  link_graph_dis;

class link_graph_wrapper : public StoreGraph {
 private:
  link_graph_dis lgd;
  GPartitioner g_part;
 public:
  /// Adds an edge to the graph
  void add(int i, int j) { lgd.add(i,j); }
  /// callback function, returns the set of neighbors to #j# vertex. 
  void set_ngbrs(int j,GSet &ngbrs_v) { lgd.set_ngbrs(j,ngbrs_v); }
  /// Clean all memory related 
  ~link_graph_wrapper() { lgd.clear(); };
  /// Constructor
  link_graph_wrapper(int N=0,const DofPartitioner *pp=NULL,
		     MPI_Comm comm_=MPI_COMM_WORLD) :
    g_part(pp),
    lgd(&g_part,comm_) { init(N); }
  /// perform the scatter of elements to its corresponding processor. 
  void scatter() { lgd.scatter(); }
  void print() const { lgd.print(); }
  /// Clean all memory related 
  void clear() { lgd.clear(); }
};

#endif
