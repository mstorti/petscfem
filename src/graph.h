// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: graph.h,v 1.12 2002/07/18 02:57:30 mstorti Exp $
#ifndef GRAPH_H
#define GRAPH_H

#include <queue>
#include <set>

/// A set of neighbors. 
typedef set<int,less<int>,malloc_alloc> GSet;

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
  void clear();
  /// Print the graph
  void print();
};

#endif
