/*__INSERT_LICENSE__*/
// $Id: tryme3.cpp,v 1.2 2002/07/23 12:27:54 mstorti Exp $

#include <src/utils.h>
#include <src/linkgraph.h>

void row_print(LinkGraphRow row) {
  LinkGraphRow::iterator q;
  for (q=row.begin(); q!=row.end(); q++) printf("%d ",*q);
  printf("\n");
}

int main(int argc, char **args) {
  MPI_Init(&argc,&args);
  LinkGraphDis graph;
  int M=10;
  graph.init(M);
  for (int j=0; j<M; j++) {
    graph.add(j,modulo(j+1,M));
    graph.add(j,modulo(j-1,M));
  }

  LinkGraphDis::iterator k;
  for (k=graph.begin(); k!=graph.end(); k++) {
    row_print(*k);
  }
}
