/*__INSERT_LICENSE__*/
// $Id: tryme3.cpp,v 1.14 2002/07/24 17:24:52 mstorti Exp $
#define _GNU_SOURCE

#include <src/utils.h>
#include <src/linkgraph.h>

int MY_RANK,SIZE;
const int M = 30;
const int N = 5;

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
/// Size of packed row (plus header)
int LinkGraphDis::size_of_pack(Row const & row) const {
  int n = row.size();
  // size + row number + size*(int+double)
  return (n+2)*sizeof(int);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Pack the row
void LinkGraphDis::pack(const Row & row,char *&buff) const {
  int n=row.size();
  BUFFER_PACK<int>(row.row,buff);
  BUFFER_PACK<int>(n,buff);
  Row::iterator q;
  for (q=row.begin(); q!=row.end(); q++)
    BUFFER_PACK<int>(*q,buff);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinkGraphDis::unpack(Row & row,const char *&buff)  {
  int n,k;
  BUFFER_UNPACK<int>(row.row,buff);
  BUFFER_UNPACK<int>(n,buff);
  for (int j=0; j<n; j++) {
    BUFFER_UNPACK<int>(k,buff);
    row.insert(k);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// combine a row in the container
void LinkGraphDis::combine(const Row &row) {
  int j=row.row;
  Row::iterator q;
  for (q=row.begin(); q!=row.end(); q++) list_insert(j,*q);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void graph_print(LinkGraphDis &graph, char *s=NULL) {
  if (s) printf("%s\n",s);
  for (LinkGraph::iterator k=graph.begin(); k!=graph.end(); k++) {
    row_print(*k);
    // printf("size of row: %d\n",k.size());
  }
}

void graph_print_dis(LinkGraphDis &graph, char *s=NULL) {
  if (MY_RANK==0) printf("%s",s);
  for (int p=0; p<SIZE; p++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (MY_RANK==0) printf("On [%d]\n",p);
    if (p==MY_RANK) graph_print(graph);
    fflush(stdout);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main(int argc, char **args) {

  MPI_Init(&argc,&args);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);
  MPI_Comm_rank (MPI_COMM_WORLD, &MY_RANK);

#define BUFSIZE 10000
  char buff[BUFSIZE], *buffd;
  const char *buffc;

#if 0
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // SEQUENTIAL DEBUG
  LinkGraphDis graph(&part,MPI_COMM_WORLD,1),
    graph2(&part,MPI_COMM_WORLD),
    *from, *to;
  
  graph.init(M);
  graph2.init(M);
  for (int j=0; j<int(sqrt(M))*M; j++) {
    graph.add(irand(M),irand(M));
    graph2.add(irand(M),irand(M));
  }

  printf("-------------\nBefore copying: \n");
  graph_print(graph,"graph: ");
  graph_print(graph2,"graph2: ");
  
  int indx=0;
  LinkGraph::iterator q=graph.begin(), q2=graph2.begin(), qq;
  while (q!=graph.end() && q2!=graph2.end()) {
    if (indx % 2 != 0) { from = &graph; to = &graph2; qq=q; }
    else { from = &graph2; to = &graph; qq=q2; }
    LinkGraphRow row = *qq,roww;
    int n = from->size_of_pack(row);
    assert(n<=BUFSIZE);
    buffd=buff;
    from->pack(row,buffd);
    from->erase(qq);

    buffc=buff;
    to->unpack(roww,buffc);
    to->combine(roww);
    indx++; q++; q2++;
  }

  printf("-------------\nAfter copying: \n");
  graph_print(graph,"graph: ");
  graph_print(graph2,"graph2: ");

#else

  // ================================================================
  // TRY SCATTER
  LinkGraphDis graph(&part,MPI_COMM_WORLD,1);
  
  graph.init(M);
  for (int j=0; j<M; j++) {
    for (int k=-N; k<=N; k++) 
      if (modulo(j+k,SIZE)==MY_RANK || modulo(j+k+1,SIZE)==MY_RANK)
	graph.add(j,modulo(j+k,M));
  }

  graph_print_dis(graph,"-----------\nBefore scatter:\n");
  graph.scatter();
  graph_print_dis(graph,"-----------\nAfter scatter:\n");

  MPI_Finalize();
#endif

}
