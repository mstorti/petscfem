// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: graph.h,v 1.1 2001/11/22 23:07:51 mstorti Exp $
#ifndef GRAPH_H
#define GRAPH_H

class Graph {
 public:
  ~Graph()=0;
  /// The number of vertices in the fine `graph'
  int nelemfat;
  /// The number of vertices in the coarse `graph'
  int nvrtx;
  /// max number of vertices in the allowed in the coarse `graph'
  int max_partgraph_vertices;
  /// number of subdomains to partition
  int npart;
  /** callback function for the user to return the neighbors of
      fine vertex elem
  */
  void ngbrs(int elem);
  /** callback function for the user to return a neighbor to the 
      coarse graph vertex
  */
  void set_ngbr(int ngbr);
}

#endif
