//__INSERT_LICENSE__
//$Id: iisdgraph.cpp,v 1.1.2.1 2001/12/13 21:44:55 mstorti Exp $

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;

//  #define DEBUG_IISD
//  #define DEBUG_IISD_DONT_SET_VALUES

#include <map>
#include <typeinfo>

#include <src/graph.h>
#include <src/distcont.h>

#if 0
template <typename IntPartitioner,typename ValueType>
class Partitioner {
private:
  const int n;
public:
  void processor(const ValueType &p,int &nproc,int *plist) const {};
  //// code more stuff here .........
}
#endif

template <typename ImgValueType>
class IntPartitioner {
private:
  typedef pair<int,ImgValueType> ValueType;
public:
  virtual int processor(int j) const =0;
  void processor(const ValueType &p,int &nproc,int *plist) const {
    nproc = 1;
    plist[0] = processor(p.first);
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This `Graph' class has internal storage, which you can fill with
    `set'. Furthermore, it is distributed.
*/
class StoreGraph : public Graph {
private:
  /// The MPI compunicator
  MPI_Comm comm;
  /// The storage area type
  typedef map<int, set<int> > GMap;
  /// An individual set of the storage map. 
  typedef pair<int, set<int> > GRow;
  /** This is the storage area. The Graph is stored as a
      correspondence between an integer (one vertex) and a set of
      integers (those vertices that are connected to him).
  */
  typedef IntPartitioner< set<int> > Partitioner;
  DistCont <GMap,GRow,Partitioner> lgraph;
public:
  /// Clean all memory related 
  ~StoreGraph() { map.clear(); };
  /// Constructor
  StoreGraph(int N=0,IntPartitioner *pp=NULL,
	    MPI_Comm comm_=MPI_COMM_WORLD) :
  lgraph(pp,comm) { init(N); }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:   
void IISDGraph::set_ngbrs(int loc1,set<int> &ngbrs_v) {
  int pos,loc2,dof2;
  pos = loc2dof[loc1]+k1;
  while (1) {
    nodep = (Node *)da_ref(da,pos);
    if (nodep->next==-1) break;
    // loc2:= number of global dof connected to `'
    dof2 = nodep->val;
    if (k1<=dof2 && dof2<=k2 && !flag[dof2] ) {
      loc2 = dof2loc[dof2-k1];
      ngbrs_v.insert(loc2);
    }
    pos = nodep->next;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double IISDGraph::weight(int elem) {
  return 1.;
}

