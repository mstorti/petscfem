// -*- mode: C++ -*- 
//__INSERT_LICENSE__
//$Id: iisdgraph.h,v 1.3 2002/07/18 02:57:30 mstorti Exp $
#ifndef IISDGRAPH_H
#define IISDGRAPH_H

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;

//  #define DEBUG_IISD
//  #define DEBUG_IISD_DONT_SET_VALUES

#include <map>
#include <typeinfo>

#include <src/part.h>
#include <src/graph.h>
#include <src/distcont.h>

/// The storage area type
typedef map<int, GSet, less<int>, malloc_alloc > GMap;
/// An individual set of the storage map. 
typedef pair<int, GSet > GRow;
/// Partitioner for the scatter operation. 
typedef Partitioner< GSet > GPartitioner;
/// The basic container for the distributed graph class.
typedef DistCont<GMap,GRow,GPartitioner> DGMap;
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This `Graph' class has internal storage, which you can fill with
    `set'. Furthermore, it is distributed.
*/
class StoreGraph : public Graph {
 public:
 private:
  /// The MPI compunicator
  MPI_Comm comm;
  /** This is the storage area. The Graph is stored as a
      correspondence between an integer (one vertex) and a set of
      integers (those vertices that are connected to him).
  */
  DGMap lgraph;
  GPartitioner g_part;
 public:
  /// Adds an edge to the graph
  void add(int i, int j) { lgraph[i].insert(j); }
  /// callback function, returns the set of neighbors to #j# vertex. 
  void set_ngbrs(int j,GSet &ngbrs_v);
  /// Clean all memory related 
  ~StoreGraph() { lgraph.clear(); };
  /// Constructor
  StoreGraph(int N=0,const DofPartitioner *pp=NULL,
	     MPI_Comm comm_=MPI_COMM_WORLD) :
    g_part(pp),
    lgraph(&g_part,comm_), comm(comm_) { init(N); }
  // void print() { lgraph.print(); }
  /// perform the scatter of elements to its corresponding processor. 
  void scatter() { lgraph.scatter(); }
  void print() const;
};

#endif
