/*__INSERT_LICENSE__*/
// $Id: tryme.cpp,v 1.11 2002/01/14 03:45:06 mstorti Exp $

#include <src/utils.h>
#include <src/graph.h>

class TGraph : public Graph {
public:
  int N,M;
  void n2ij(int n,int &j,int &k);
  int ij2n(int j,int k);
  void set_ngbrs(int elem,set<int> &ngbrs_v);
  ~TGraph() {clear();}
};

void TGraph::n2ij(int n,int &j,int &k) {
  k = n/N;
  j = n % N;
}

int TGraph::ij2n(int j,int k) {
  return k*N+j;
}

void TGraph::set_ngbrs(int elem,set<int> &ngbrs_v) {
  int n,j,k;
  n = elem;
  n2ij(n,j,k);
  ngbrs_v.insert(ij2n(crem(j+1,N),k));
  ngbrs_v.insert(ij2n(crem(j-1,N),k));
  ngbrs_v.insert(ij2n(j,crem(k+1,M)));
  ngbrs_v.insert(ij2n(j,crem(k-1,M)));

#if 0
  printf("node %d adding %d %d %d %d\n",n,ij2n(crem(j+1,N),k),
	 ij2n(crem(j-1,N),k), ij2n(j,crem(k+1,M)),
	 ij2n(j,crem(k-1,M)));
#endif
  
}

int main(int argc, char **args) {
  int N=200,M=1;
  TGraph GG;
  Graph &G = GG;
  int j,k=0,compact=0,arg=0,npart=2;

  if (argc>arg++) {
    sscanf(args[arg],"%d",&N);
  }
  assert(N>0);

  if (argc>arg++) {
    sscanf(args[arg],"%d",&M);
  }
  assert(M>0);

  if (argc>arg++) {
    sscanf(args[arg],"%d",&npart);
  }
  assert(M>0);

  G.init(N*M);

  if (argc>arg++) {
    sscanf(args[arg],"%d",&compact);
  }

  if (compact) {
    assert(M<=60);
    assert(N<=30);
    assert(npart<=10);
  }

#if 1
  printf("partitioning graph %d\n",k++);
  GG.N = N;
  GG.M = M;
  G.part(N*M,npart);
  if (compact) {
    for (j=0; j<N; j++) {
      for (k=0; k<M; k++)
	printf("%1d",G.vrtx_part(GG.ij2n(j,k)));
      printf("\n");
    }
  } else {
    for (int j=0; j<N; j++) 
      for (int k=0; k<M; k++) 
	printf("vrtx %d %d in proc %d\n",
	       j,k,G.vrtx_part(GG.ij2n(j,k)));
  }
  G.clear();
#endif

#if 0 // debugging memory leakages
  while (1) {
    printf("partitioning graph %d\n",k++);
    GG.N = N;
    G.part(N,N,2);
    G.clear();
//      for (int j=0; j<N; j++) 
//        printf("vrtx %d in proc %d\n",j,G.vrtx_part(j));
  }
#endif
}
