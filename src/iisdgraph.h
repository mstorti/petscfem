//__INSERT_LICENSE__
//$Id: iisdgraph.h,v 1.1.2.1 2001/12/17 00:03:57 mstorti Exp $

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;

//  #define DEBUG_IISD
//  #define DEBUG_IISD_DONT_SET_VALUES

#include <map>
#include <typeinfo>

#include <src/part.h>
#include <src/graph.h>
#include <src/distcont.h>

template <typename ImgValueType>
class Partitioner  {
private:
  typedef pair<int,ImgValueType> ValueType;
public:
  DofPartitioner *part;
  void processor(const ValueType &p,int &nproc,int *plist) const {
    nproc = 1;
    plist[0] = part->processor(p.first);
  }
  Partitioner<ImgValueType> (DofPartitioner *pp) 
    : part(pp) { assert(pp!=NULL) ; }
};

/// The storage area type
typedef map<int, set<int> > GMap;
/// An individual set of the storage map. 
typedef pair<int, set<int> > GRow;
/// Partitioner for the scatter operation. 
typedef Partitioner< set<int> > GPartitioner;
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
 public:
  /// Adds an edge to the graph
  void add(int i, int j) { lgraph[i].insert(j); }
  /// callback function, returns the set of neighbors to #j# vertex. 
  void set_ngbrs(int j,set<int> &ngbrs_v);
  /// Clean all memory related 
  ~StoreGraph() { lgraph.clear(); };
  /// Constructor
  StoreGraph(int N=0,GPartitioner *pp=NULL,
	     MPI_Comm comm_=MPI_COMM_WORLD) :
    lgraph(pp,comm_) { init(N); }
  // void print() { lgraph.print(); }
  /// perform the scatter of elements to its corresponding processor. 
  void scatter() { lgraph.scatter(); }
#if 0
  /** Computes the size of data needed to pack this entry 
      @param k (input) iterator to the entry
      @return the size in bytes of the packed object
  */ 
  int size_of_pack(const GRow &p) const;
  /** Packs the entry #(k,v)# in buffer #buff#. This function should
      be defined by the user. 
      @param k (input) key of the entry
      @param v (input) value of the entry
      @param buff (input/output) the position in the buffer where the
      packing is performed
  */ 
  void pack(const GRow &p,char *&buff) const;
  /** Does the reverse of #pack#. Given a buffer #buff# recovers the
      corresponding key and val. This function should
      be defined by the user. 
      @param k (output) key of the entry
      @param v (output) value of the entry
      @param buff (input/output) the position in the buffer from where the
      unpacking is performed
  */ 
  void unpack(GRow &p,const char *& buff);
#endif
};

