/*__INSERT_LICENSE__*/
// $Id: tryme3.cpp,v 1.5 2002/07/24 01:19:58 mstorti Exp $

#include <src/utils.h>
#include <src/linkgraph.h>

int MY_RANK,SIZE;
const int M=10;

void row_print(LinkGraphRow row) {
  LinkGraphRow::iterator q;
  printf("row %d: ",row.row);
  for (q=row.begin(); q!=row.end(); q++) printf("%d ",*q);
  printf("\n");
}

class Part {
public:
  int processor(int j) { return int((j*SIZE)/M);};
  void processor(const LinkGraphRow &k,int &nproc,int *plist);
} part;

void 
Part::processor(const LinkGraphRow &k,int &nproc,int *plist) {
  nproc=1;
  plist[0] = processor(k.row);
}

typedef  
DistCont<LinkGraph,LinkGraphRow,Part>  LinkGraphDis;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main(int argc, char **args) {
  MPI_Init(&argc,&args);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank (MPI_COMM_WORLD, &MY_RANK);

  LinkGraphDis graph(&part,MPI_COMM_WORLD);
  graph.init(M);
  for (int j=0; j<M; j++) {
    graph.add(j,modulo(j+1,M));
    graph.add(j,modulo(j-1,M));
  }
  graph.scatter();
  LinkGraphDis::iterator k;
  for (k=graph.begin(); k!=graph.end(); k++) {
    row_print(*k);
    printf("size of row: %d\n",k.size());
  }
}
