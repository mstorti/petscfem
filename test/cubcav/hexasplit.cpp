/*__INSERT_LICENSE__*/
// $Id: hexasplit.cpp,v 1.2 2002/07/28 20:51:36 mstorti Exp $
#define _GNU_SOURCE

#include <src/utils.h>
#include <src/linkgraph.h>
#include <src/dvector.h>

int MY_RANK,SIZE;

void row_print(LinkGraphRow row) {
  LinkGraphRow::iterator q;
  printf("row %d: ",row.row);
  for (q=row.begin(); q!=row.end(); q++) printf("%d ",*q);
  printf("\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void graph_print(LinkGraphDis &graph, char *s=NULL) {
  if (s) printf("%s\n",s);
  for (LinkGraph::iterator k=graph.begin(); k!=graph.end(); k++) {
    row_print(*k);
    // printf("size of row: %d\n",k.size());
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main(int argc, char **args) {

  // SEQUENTIAL DEBUG
  LinkGraph graph;
  // graph.init(M);
  char *line=NULL;
  size_t ll=0;
#define NEL 8
  int nodes[NEL],node,nnod=0;
  dvector<int> icone;
#define MIN_CHUNK_SIZE 40000
  icone.set_chunk_size(MIN_CHUNK_SIZE);
  // Incompatible nodes in a hexa
  int incompat[] = {0,1,0,1,1,0,1,0};

  FILE *fid = fopen(args[1],"r");
  assert(fid);
  int nelem=0;
  while(1) {
    char *token;
    int linelen = getline(&line,&ll,fid);
    if (linelen<0) break;
    nelem++;
    for (int k=0; k<NEL; k++) {
      token = strtok((k==0? line : NULL)," ");
      int nread = sscanf(token,"%d",&node);
      assert(nread==1);
      icone.push(node);
      if (node>nnod) nnod=node;
    }
  }
  printf("read %d elems, %d nodes\n",nelem,nnod);

  fclose(fid);
  graph.set_chunk_size(nnod/2 < MIN_CHUNK_SIZE ? MIN_CHUNK_SIZE : nnod/2);
  graph.init(nnod);
  for (int e=0; e<nelem; e++) {
    int *row = &icone.ref(e*NEL);
    for (int j=0; j<NEL; j++) 
      for (int k=0; k<NEL; k++) 
	if (incompat[j]!=incompat[k]) graph.add(row[j]-1,row[k]-1);
  }
  for (int q=0; q<nnod; q++) {
    GSet row;
    graph.set_ngbrs(q,row);
    printf("row %d: ",q);
    for (GSet::iterator r=row.begin(); r!=row.end(); r++) 
      printf("%d ",*r);
    printf("\n");
  }
}
