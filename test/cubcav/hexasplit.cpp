/*__INSERT_LICENSE__*/
// $Id: hexasplit.cpp,v 1.1 2002/07/28 20:30:17 mstorti Exp $
#define _GNU_SOURCE

#include <src/utils.h>
#include <src/linkgraph.h>

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
      if (node>nnod) nnod=node;
      assert(nread==1);
    }
  }
  printf("read %d elems, %d nodes\n",nelem,nnod);
}
