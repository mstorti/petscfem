// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: graph.h,v 1.3 2001/11/23 02:13:28 mstorti Exp $
#ifndef GRAPH_H
#define GRAPH_H

#include <queue>
#include <vector>
#include <set>

class Graph {
  /*** ngbrs:= Auxiliary queue that stores those nodes that are waiting
       to be marked. 
  */
  queue<int> ngbrs; 
  // adjncy:= xadj:= graph described in CSR format (as defined in
  // Metis documentation)
  int *xadj;
  // vwgt:= weights for the vertices (elements) of the graph
  int *vwgt;
  // el2vrtx:= maps elements to vertices when coalescing
  int *el2vrtx;
  // vpart:= partitioning of vertices
  int *vpart;
  // Stores the graph `adjncy' in STL format
  vector< set<int> > adjncy_v;
  set<int> subd_remap_aux;
  double weight_scale;
  /// The number of vertices in the fine `graph'
  int nelemfat;
  /// The number of vertices in the coarse `graph'
  int nvrtx;
  /// max number of vertices in the allowed in the coarse `graph'
  int max_partgraph_vertices;
 public:
  /// Constructor
  Graph(int nelemfat,int max_partgraph_vertices);
  /** Virtual destructor (Just to prevent instantiations of this
      class. See Thinking in C++, by Bruce Eckel).
  */
  virtual ~Graph()=0;
  /// number of subdomains to partition
  void part(int npart);
  /// return an array with the domain indices of all fine vertices 
  const int *vrtx_part();
  /// return the partition index for a given fine graph vertex 
  int vrtx_part(int);
  /** callback user function to return the neighbors for a 
      fine vertex `elem'
  */
  virtual void set_ngbrs(int elem,vector<int> &ngbrs_v);
  /// Callback user function for the user to set the weight of a given fine vertex. 
  double weight(int elem);
  /// Clean all memory related 
  void clear();
};

#endif
