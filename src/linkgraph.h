// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: linkgraph.h,v 1.8 2002/07/23 17:06:59 mstorti Exp $
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
  void erase(iterator q) { assert(0); } // not coded yet...
  void erase(iterator first,iterator last) { assert(0); } // not coded yet...
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
#if 0
  int size_of_pack(const_iterator iter) const {
    int n = iter->second.size();
    // size + row number + size*(int+double)
    return (n+2)*sizeof(int)+n*sizeof(double);
  }
#endif
};

typedef LinkGraph::Row LinkGraphRow; 
// typedef Partitioner<LinkGraphRow> LinkGraphRowPart;

//  typedef  
//  DistCont<LinkGraph,LinkGraphRow,LinkGraphRowPart>  LinkGraphDis;

#if 0
class LinkGraphWrapper : public StoreGraph {
 private:
  LinkGraphDis lgd;
  GPartitioner g_part;
 public:
  /// Adds an edge to the graph
  void add(int i, int j) { lgd.add(i,j); }
  /// callback function, returns the set of neighbors to #j# vertex. 
  void set_ngbrs(int j,GSet &ngbrs_v) { lgd.set_ngbrs(j,ngbrs_v); }
  /// Clean all memory related 
  ~LinkGraphWrapper() { lgd.clear(); };
  /// Constructor
  LinkGraphWrapper(int N=0,const DofPartitioner *pp=NULL,
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

#endif
