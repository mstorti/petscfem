// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: linkgraph.h,v 1.14 2003/07/02 23:22:19 mstorti Exp $
#ifndef LINKGRAPH_H
#define LINKGRAPH_H

extern int MY_RANK,SIZE;

#include <src/buffpack.h>
#include <src/distcont.h>
#include <src/graph.h>
#include <src/iisdgraph.h>
#include <src/dvector.h>
#include <src/graphdv.h> // for int_pair

/** Graph adjacency class that stores the list of vertices 
    adjacent to a vertex #i# as a linked list of cells pointed 
    by cursors in a dynamic vector (dvector). 
*/
class LinkGraph {
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
  /// free a cell (put it in free list)
  void free_cell(int cell);
public:
  class Row : public GSet {
  public:
    int row;
  };
  // class Row : public set<int> {};

  class iterator;
  friend class iterator;
  class const_iterator {
    friend class LinkGraph;
  protected:
    int r;
    LinkGraph *graph;
    Row row_;
  public:
    const_iterator(int rr=0,LinkGraph *g=NULL) : r(rr), graph(g) {}
    const_iterator &operator++(int) { 
      r++;
      return *this;
    }
    const Row& operator*() {
      row_.row = r;
      graph->set_ngbrs(r,row_);
      return row_;
    }
    int operator==(const_iterator q) const { return r==q.r; }
    int operator!=(const_iterator q) const { return r!=q.r; }
    /// returns number of adjacents elements in the row
    int size() const { graph->size(r); }
  };
  class iterator : public const_iterator {
    friend class LinkGraph;
  public:
    iterator(int rr=0,LinkGraph *g=NULL) : const_iterator(rr,g) {}
    Row operator*() {
      Row row;
      row.row = r;
      graph->set_ngbrs(r,row);
      return row;
    }
  };    
  iterator begin() { return iterator(0,this); }
  iterator end() { return iterator(M,this); }
  void erase(iterator q);
  void erase(iterator first,iterator last);
  /// cuasi constructor
  void init(int MM);
  /** Set new chunk size (container must be empty).
      @param new_chunk_size (input) new chunk size for the vector.
  */
  void set_chunk_size(int new_chunk_size) { da.set_chunk_size(new_chunk_size); }
  /// Destructor. Clear dynamic array. 
  ~LinkGraph() { da.clear(); }
  /** Add an edge to the adjacency table
    @param j (input) first index of edge
    @param k (input) second index of edge
  */
  void add(int j,int k) { list_insert(j,k); }
  /** Queries whether a given pair of indices
      is linked by an edge or not. 
      @param j (input) first vertex
      @param k (input) second vertex
      @return boolean value indicating if the pair is an edge or not. 
  */ 
  int edge_q(int j,int k);
  /** Call back function that defines the graph. For a given
      vertex it returns the list of adjacent vertices. 
      @param v (input) the source vertex 
      @param ngbrs (output) the set of adjacent vertices 
  */
  void set_ngbrs(int v,GSet &ngbrs);
  /// Return number of edges
  int size() { return da.size(); } // not implemented yet
  /// number of vertices adjacent to vertex r
  int size(int r);
  /// Clear all edges
  void clear() { M=0; da.clear(); }
};

typedef LinkGraph::Row LinkGraphRow; 

//typedef Partitioner<LinkGraphRow> LinkGraphRowPart;
class LinkGraphRowPart {
private:
  const DofPartitioner *part;
public:
  LinkGraphRowPart(const DofPartitioner *dp) : part(dp) { assert(dp); }
  void processor(const LinkGraphRow &p,int &nproc,int *plist) const {
    nproc = 1;
    plist[0] = part->processor(p.row);
  }
};

typedef  
DistCont<LinkGraph,LinkGraphRow,LinkGraphRowPart>  LinkGraphDis;
class LinkGraphWrapper : public StoreGraph {
 private:
  LinkGraphDis lgd;
  LinkGraphRowPart lg_part;
 public:
  /// Adds an edge to the graph
  void add(int i, int j) { lgd.add(i,j); }
  /// callback function, returns the set of neighbors to #j# vertex. 
  void set_ngbrs(int j,GSet &ngbrs_v) { lgd.set_ngbrs(j,ngbrs_v); }
  /// Clean all memory related 
  ~LinkGraphWrapper() { lgd.clear(); };
  /// Constructor
  LinkGraphWrapper(int N=0,const DofPartitioner *dp=NULL,
		   MPI_Comm comm_a=MPI_COMM_WORLD) :
    lg_part(dp),
    lgd(&lg_part,comm_a,LinkGraphDis::random_iter_mode) { init(N); }
  /// perform the scatter of elements to its corresponding processor. 
  void scatter() { lgd.scatter(); }
  /// Clean all memory related 
  void clear() { lgd.clear(); }
  /** Set new chunk size.
      @param new_chunk_size (input) new chunk size for the vector.
  */
  void set_chunk_size(int new_chunk_size) { lgd.set_chunk_size(new_chunk_size); }
  /// 
  void init(int MM) { Graph::init(MM); lgd.init(MM); }
};

#endif
