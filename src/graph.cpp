//__INSERT_LICENSE__
//$Id: graph.cpp,v 1.8 2001/11/26 22:41:43 mstorti Exp $

#include <src/utils.h>
#include <src/graph.h>
extern "C" {
#include <metis.h>
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Graph::~Graph() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Graph::Graph() : weight_scale(1.), vpartf(NULL) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double Graph::weight(int vrtx_f) {return weight_scale;}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Graph::part(int nvrtx_f,int max_partgraph_vertices,
  int npart,float *tpwgts=NULL) {

  int visited,vrtx_f,vrtx_fk,vrtx,vrtxj,vrtxjj,p,j,k,
    edgecut,options=0,numflag=0,wgtflag=2;
  vector<int>::iterator n,ne;
  set<int>::iterator q,qe;
  // if tpwgts is not passed then define a local one
  vector<float> tpwgts_v;
  // Stores the graph `adjncy' in STL format
  vector< set<int> > adjncy_v;

  weight_scale = 1.;

  assert(nvrtx_f>0);
  assert(max_partgraph_vertices>0);
  
  // nvrtx:= number of vertices used in graph partitioning
  // The factor 2 here avoids the case that only few elements are not
  // coalesced. In this case, for the last elements the computational
  // effort may be high, since the probability of getting a hole is
  // very low. 
  nvrtx = (nvrtx_f/2 > max_partgraph_vertices ?
	   max_partgraph_vertices : nvrtx_f);
  
  // Create adjacency table for partitioning with Metis. In the
  // adjacency graph the nodes are elements or group of elements of the
  // FEM mesh. Two nodes of the graph (elements of the mesh) have are
  // linked if they share a node.

  // adjncy:= xadj:= graph described in CSR format (as defined in
  // Metis documentation)
  int *xadj = new int[nvrtx+1],*adjncy;
  // vwgt:= weights for the vertices (elements) of the graph
  int *vwgt = new int[nvrtx];
  // el2vrtx:= maps elements to vertices when coalescing
  el2vrtx = new int[nvrtx_f];
  // vpart:= partitioning of vertices
  vpart = new int[nvrtx];
  vector<int> ngbrs_v;

  // Initialize 
  for (j=0; j<nvrtx_f; j++) el2vrtx[j] = -1;
  for (j=0; j<nvrtx; j++) vwgt[j] = 0;
  // Take `nvrtx' elements to coalesce the neighbor elements

  vrtx=0;
  while (!ngbrs.empty()) ngbrs.pop();

  if (nvrtx_f!=nvrtx) {
    // Take nvrtx `seeds' for growing the groups 
    while(1) {
      vrtx_f = irand(nvrtx_f);
      if (el2vrtx[vrtx_f]!=-1) continue;
      el2vrtx[vrtx_f]=vrtx++;
      vwgt[vrtx] += int(weight(vrtx_f)/weight_scale);
      ngbrs.push(vrtx_f);
      if (vrtx==nvrtx) break;
    }
    visited = nvrtx;
  } else {
    // Take a group for each element
    for (vrtx_f=0; vrtx_f<nvrtx_f; vrtx_f++) {
      el2vrtx[vrtx_f] = vrtx_f;
      vwgt[vrtx_f] += int(weight(vrtx_f)/weight_scale);
    }      
    visited=nvrtx_f;
  }
  
  // Start coalescing neighbor elements until all elements are
  // assigned a vertex
  while (1) {
    if (ngbrs.empty()) {
      // either we have visited all the elements or there are
      // disconnected parts of the domain.

      // All elements have been visited
      if (visited==nvrtx_f) break;
      
      // Disconnected parts. Find next element not visited.
      for (k=0; k<nvrtx_f; k++) {
	if (el2vrtx[k]<0) {
	  ngbrs.push(k);
	  // Put it in a random vertex
	  vrtx = irand(0,nvrtx)-1;
	  el2vrtx[k] = vrtx;
	  vwgt[vrtx] += int(weight(vrtx_f)/weight_scale);
	  break;
	}
      }
    }

    // Take the first element
    vrtx_f = ngbrs.front();
    ngbrs.pop();
    vrtx = el2vrtx[vrtx_f];
    assert(vrtx>=0);
    // Get its connectivities
    ngbrs_v.clear();
    set_ngbrs(vrtx_f,ngbrs_v);
    // Loop over all ngbrs of the fine vertex
    ne = ngbrs_v.end();
    for (n=ngbrs_v.begin(); n!=ne; n++) {
      int &vrtx_fk = *n;
      // Already visited
      if (el2vrtx[vrtx_fk]>=0) continue;
      // Assign to group
      el2vrtx[vrtx_fk]=vrtx;
      visited++;
      // Find elemset and add weight to vertex weight (for load
      // balancing) 
      vwgt[vrtx] += int(weight(vrtx_f)/weight_scale);
      // Put in queue
      ngbrs.push(vrtx_fk);
    }
  }

#if 0
  if (myrank==0) {
    for (vrtx_f=0; vrtx_f < nvrtx_f; vrtx_f++) {
      printf("el2vrtx[%d] = %d\n",vrtx_f,el2vrtx[vrtx_f]);
    }
  }
#endif

  // Resize the `vector of sets'
  adjncy_v.resize(nvrtx);

  // Define graph descriptors 
  // Initializa graph ptrs
  for (vrtx_f=0; vrtx_f<nvrtx_f; vrtx_f++) {
    vrtxj = el2vrtx[vrtx_f];
    ngbrs_v.clear();
    set_ngbrs(vrtx_f,ngbrs_v);
    ne = ngbrs_v.end();
    for (n=ngbrs_v.begin(); n!=ne; n++) {
      int &vrtx_fk = *n;
      vrtxjj = el2vrtx[vrtx_fk];
      // Add edges to graph
      adjncy_v[vrtxj].insert(vrtxjj);
      // Just in case the user doesn't return a symmetric graph
      adjncy_v[vrtxjj].insert(vrtxj);
    }
  }

  xadj[0] = 0;
  for (vrtxj=0; vrtxj<nvrtx; vrtxj++) 
    xadj[vrtxj+1] = xadj[vrtxj] + adjncy_v[vrtxj].size();

  // Once computed the adjncy size, it is created.
  adjncy = new int[xadj[nvrtx]];
  for (vrtxj=0; vrtxj<nvrtx; vrtxj++) {
    set<int> &adj = adjncy_v[vrtxj];
    qe = adj.end();
    p = xadj[vrtxj];
    for (q=adj.begin(); q!=qe; q++) adjncy[p++] = *q;
  }

  adjncy_v.clear();
#if 0
  // print the graph
  for (int ielgj=0; ielgj<nvrtx; ielgj++) {
    printf("%d: ",ielgj);
    for (int jj=xadj[ielgj]; jj<xadj[ielgj+1]; jj++) 
      printf(" %d",adjncy[jj]);
    printf(", weight: %d\n",vwgt[ielgj]);
  }
#endif

  if (tpwgts==NULL) {
    tpwgts_v.resize(npart);
    for (k=0; k<npart; k++) tpwgts_v[k] = 1./float(npart);
    tpwgts = tpwgts_v.begin();
  }

  if (npart>1) {
    METIS_WPartGraphKway(&nvrtx,xadj,adjncy,vwgt, 
			 NULL,&wgtflag,&numflag,&npart, 
			 tpwgts,&options,&edgecut,vpart);
  } else {
    for (j=0; j<nvrtx; j++) vpart[j] = 0;
  }

  delete[] adjncy;
  delete[] xadj;
  delete[] vwgt;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Graph::vrtx_part(int vrtx_f) {
  // Map fine graph vertex to coarse graph vertex and return partition
  // index for the last one
  return vpart[el2vrtx[vrtx_f]];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
const int *Graph::vrtx_part() {
  // Create auxiliary vector `vpartf', fill with values using
  // elemental `vrtx_part' and return pointer to it
  vpartf = new int[nvrtx_f];
  for (int j=0; j<nvrtx_f; j++) 
    vpartf[j] = vpart[el2vrtx[j]];
  return vpartf;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void Graph::clear() {
  delete[] el2vrtx;
  el2vrtx = NULL;
  delete[] vpart;
  vpart = NULL;
  delete[] vpartf;
  vpartf = NULL;
}
