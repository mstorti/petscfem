/*__INSERT_LICENSE__*/
// $Id: pfmat.cpp,v 1.1 2001/12/08 20:31:00 mstorti Exp $

// Tests for the `PFMat' class
#include <src/utils.h>
#include <src/util2.h>
#include <src/graph.h>
#include <src/pfmat.h>

static char help[] = "PETSc-FEM Navier Stokes module\n\n";

/// A graph that has an internal representation
class IntGraph : public Graph {
private:
  static const int END = -1;
  Darray *da;
public:
  IntGraph(int N=0) { init(N); };
  void init(int N);
  ~IntGraph();
  void set_ngbrs(int vrtx_f,set<int> &ngbrs_v);
  void add_ngbr(int v1,int v2);
};

IntGraph::~IntGraph() {}

void IntGraph::init(int N) {
  if (N==0) return;
  Graph::init(N);
  da= da_create_len(sizeof(Node),N);
  Node nodeq = Node(END,0);
  for (int k=0; k < N; k++) {
    da_set(da,k,&nodeq);
  }
}

void IntGraph::set_ngbrs(int vrtx_f,set<int> &ngbrs_v) {
  Node *nodep,nodeq;
  int posit=vrtx_f,newpos;
  while (1) {
    nodep = (Node *)da_ref(da,posit);
    if (nodep->next != END) {
      ngbrs_v.insert(nodep->val);
    } else {
      break;
    } 
    posit = nodep->next;
  }
}

void IntGraph::add_ngbr(int j,int k) {
  Node *nodep,nodeq;
  int posit=j,newpos;

  //  printf(" Inserting %d %d, pasing by (",j,k);
  while (1) {
    nodep = (Node *)da_ref(da,posit);
    if (nodep->next == END) {
      nodeq = Node(END,0);
      newpos = da_append(da,&nodeq);
      //      printf("). Appended at position %d\n",posit);
      nodeq = Node(newpos,k);
      da_set(da,posit,&nodeq);
      break;
    } else if (nodep->val==k) {
      //      printf("). Found at position %d\n",posit);
      break;
    } else {
      //      printf(" %d",nodep->val);
      posit = nodep->next;
    }
  }
}

int main(int argc,char **args) {

  PetscInitialize(&argc,&args,(char *)0,help);
  const int N=10;
  int j;
  IntGraph g(N);
  
  for (j=0; j<N; j++) {
    g.add_ngbr(j,crem(j-1,N));
    g.add_ngbr(j,crem(j+1,N));
  }
  g.print();
}
