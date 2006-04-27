// -*- mode: C++ -*- 
//__INSERT_LICENSE__
//$Id: iisdgraph.h,v 1.9.40.1 2006/04/27 20:31:10 rodrigop Exp $
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
//#include <src/graphdv.h>
extern MPI_Comm PETSC_COMM_WORLD;

/// The storage area type
typedef map<int, GSet, less<int> > GMap;
/// An individual set of the storage map. 
typedef pair<int, GSet > GRow;
/// Partitioner for the scatter operation. 
typedef Partitioner< GSet > GPartitioner;
/// The basic container for the distributed graph class.
typedef DistCont<GMap,GRow,GPartitioner> DGMap;

class StoreGraph : public Graph {
public:
  /// Adds an edge to the graph
  virtual void add(int i, int j)=0;
  virtual void scatter()=0;
};

template<>
int DGMap::size_of_pack(const GRow &q) const;

template<>
void DGMap::pack(const GRow &p,char *&buff) const;

template<>
void DGMap::unpack(GRow &p,const char *&buff);

template<>
void DGMap::combine(const GRow &p);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This `Graph' class has internal storage, which you can fill with
    `set'. Furthermore, it is distributed.
*/
class StoreGraph1 : public StoreGraph {
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
  ~StoreGraph1() { lgraph.clear(); };
  /// Constructor
  StoreGraph1(int N=0,const DofPartitioner *dp=NULL,
	     MPI_Comm comm_=PETSC_COMM_WORLD) :
    g_part(dp),
    lgraph(&g_part,comm_), comm(comm_) { init(N); }
  // void print() { lgraph.print(); }
  /// perform the scatter of elements to its corresponding processor. 
  void scatter() { lgraph.scatter(); }
  void print() const;
  /// Clean all memory related 
  void clear();
};

class StoreGraph2 : public GMap {
public:
  StoreGraph2(int N=0,const DofPartitioner *pp=NULL,
    MPI_Comm comm_=PETSC_COMM_WORLD) {}
  void add(int i,int j) { (*this)[i].insert(j); }
  void set_ngbrs(int j,GSet &ngbrs_v) {}
  ~StoreGraph2() { clear(); };
  void scatter() {}
  void print() {}
};

#endif
