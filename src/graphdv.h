// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: graphdv.h,v 1.6 2004/10/22 17:33:53 mstorti Exp $
#ifndef GRAPHDV_H
#define GRAPHDV_H

extern int MY_RANK,SIZE;

#include <src/iisdgraph.h>
#include <src/dvector.h>
#include <src/dvector2.h>

class int_pair { 
public: 
  int i,j; 
  int_pair(int ii=0,int jj=0) { i=ii; j=jj; }
  inline bool operator<(int_pair q) const {
    if (i!=q.i) return i<q.i;
    else return j<q.j;
  }
  inline bool operator>=(int_pair q) const {
    return !(q < *this);
  }
  inline bool operator==(int_pair q) const {
    return i==q.i && j==q.j;
  }
  inline bool operator>(int_pair q) const {
    if (i!=q.i) return i>q.i;
    else return j>q.j;
  }
  inline bool operator!=(int_pair q) const {
    if (i!=q.i) return 0;
    else return j!=q.j;
  }
};

/** Graph adjacency class that stores the #(i,j)# pairs in a sorted
    (lexicographically) dynamic array. It may be slower but is optimal
    in allocating memory. When adding elements they are stored in the
    back of the vector (if not present in the ordered part) and when
    they reach some threshold the vector is sorted and repeated
    elements deleted. 
    @author M. Storti
*/
class graphdv : public StoreGraph {
protected:
  /// Storage array
  dvector<int_pair> da;
  /// Size of the ordered part
  int ordered;
  /// Initial value for the #max# parameter
  int MAX_INIT;
  /// will trigger automatic sort if size passes #max#
  int max;
  /// flag whether new elements have been inserted or not
  int modif;
  /// sort array and eliminate repeated indices
  void resync();
  /// print array
  void print2dv();
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
};

/// Distributted graph.
class graphdv_dis : public graphdv  {
public:
  /// Constructor includes partitioner and communicator
  graphdv_dis(int M, const DofPartitioner *part,
	      MPI_Comm comm,int chunk_size=CHUNK_SIZE_DEF) : graphdv(chunk_size) {
    init(M);
  }
  /// Scatter among processors
  void scatter();
};

#endif
