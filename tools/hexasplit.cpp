/*__INSERT_LICENSE__*/
// $Id: hexasplit.cpp,v 1.5 2002/07/28 23:31:38 mstorti Exp $
#define _GNU_SOURCE

#include <vector>
#include <deque>

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
int main (int argc, char **argv) {
  char c;
  string icone_file = "icone";
  string icone_tetra = "icone_tetra";
  while ((c = getopt (argc, argv, "i:o:")) != -1) {
    switch (c) {
    case 'i':
      icone_file = string(optarg);
      break;
    case 'o':
      icone_tetra = string(optarg);
      break;
    default:
      abort ();
    }
  }

  LinkGraph graph;
  // graph.init(M);
  char *line=NULL;
  size_t ll=0;
#define NEL 8
  int nodes[NEL],node,nnod=0;
  dvector<int> icone;
  deque<int> front;

#define MIN_CHUNK_SIZE 40000
  icone.set_chunk_size(MIN_CHUNK_SIZE);
  // Incompatible nodes in a hexa. Two nodes have the same entry
  // (0/1) if they are not linked by an edge. 
  int incompat[] = {0,1,0,1,1,0,1,0};

  // Reads connectivities
  FILE *fid = fopen(icone_file.c_str(),"r");
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
  fclose(fid);
  printf("read %d elems, %d nodes\n",nelem,nnod);
  // split[j] may be -1/+1 depending on whether the
  // node is marked up or down. split[j]==0 implies
  // that the node is not split yet. 
  vector<int> split(nnod);
  for (int j=0; j<nelem; j++) split[j]=0;

  // Build the graph of incompatibilities. Two nodes are connected
  // by an edge if they are connected by an edge of an hexa. 
  graph.set_chunk_size(nnod/2 < MIN_CHUNK_SIZE ? MIN_CHUNK_SIZE : nnod/2);
  graph.init(nnod);
  for (int e=0; e<nelem; e++) {
    int *row = &icone.ref(e*NEL);
    for (int j=0; j<NEL; j++) 
      for (int k=0; k<NEL; k++) 
	// Nodes are incompatible if they have not the same entry in
	// the incompatibility table
	if (incompat[j]!=incompat[k]) graph.add(row[j]-1,row[k]-1);
  }

#if 0
  // Print the graph
  for (int q=0; q<nnod; q++) {
    GSet row;
    graph.set_ngbrs(q,row);
    printf("row %d: ",q);
    for (GSet::iterator r=row.begin(); r!=row.end(); r++) 
      printf("%d ",*r);
    printf("\n");
  }
#endif

  // GREEDY COLORING ALGORITHM
  // Starting from an arbitrary node we put it in a queue. Then at a
  // turn we take a node from the queue, put all their not-yet-marked
  // neighbors in the queue and mark it as up or down depending on the
  // color of their neighbors. 

  // Arbitrarily start from node `0'.
  // Split start node as `up'
#define QUEUED (-2)
  front.push_back(0);
  split[0]=1;
  while (front.size()) {
    // printf("front.size %d\n",front.size());
    // Take same node from the queue
    int node = front.front();
    front.pop_front();
    // `s' is the set of neighbors of `node'. 
    GSet s;
    graph.set_ngbrs(node,s);
    // iterate on the neighbors of `node'
    for (GSet::iterator r=s.begin(); r!=s.end(); r++) {
      if (split[*r]==0) {
	// If not already split put it in the queue
	front.push_back(*r);
	split[*r]=QUEUED;
      } else if (split[*r]==QUEUED) {
	// do nothing
      } else if (split[node]==0 || split[node]==QUEUED) {
	// if not already marked, mark as the opposite of the neighbor
	split[node] = -split[*r];
      } else if (split[*r] == split[node]) {
	// if already marked and find an incompatible neighbor complain
	printf("Can't split the mesh!!\n");
	exit(1);
      }
    }
  }

#if 0
  for (int k=0; k<nnod; k++) printf("split[%d] = %d\n",k,split[k]);
#endif

  // Print results. We simply look at the value of split[]
  // at some node of the element (say node `0') and depending of
  // this value we take the standard split on a remapped element
  // connectivity. `map_up' is the identity map, and `map_down' is
  // a map rotated 90 degrees. 
  int map_up[] = {0,1,2,3,4,5,6,7};
  int map_down[] = {1,2,3,0,5,6,7,4};
  // This will point to `map_up' or `map_down'
  int *map;
  // Standard split of an hexa in tetras. 
  int tetra[][4]={{0,2,5,1},
		  {0,7,2,3},
		  {2,7,5,6},
		  {0,5,7,4},
		  {0,5,2,7}};
  fid = fopen(icone_tetra.c_str(),"w");
  for (int k=0; k<nelem; k++) {
    // Connectivity row
    int *row = &icone.ref(k*NEL);
    // Choose map depending on split value
    map = (split[row[0]]==1 ? map_up : map_down);
    // loop over local tetras
    for (int t=0; t<5; t++) {
      // loop over nodes in the tetra
      for (int q=0; q<4; q++) 
	// node of local tetra, eventually remapped
	fprintf(fid,"%d ",row[map[tetra[t][q]]]);
      fprintf(fid,"\n");
    }
  }
  fclose(fid);
}
