//__INSERT_LICENSE__
//$Id: graph.cpp,v 1.1 2001/11/23 00:19:53 mstorti Exp $

#include <src/utils.h>
#include <src/graph.h>

void Graph::part(int npart) {
  int visited,elem,elemk,vrtx,vrtxj,vrtxjj,p;
  vector<int>::iterator n,ne;
  set<int>::iterator q,qe;

  weight_scale = 1.;

  assert(nelemfat>0);
  assert(max_partgraph_vertices>0);
  
  // nvrtx:= number of vertices used in graph partitioning
  // The factor 2 here avoids the case that only few elements are not
  // coalesced. In this case, for the last elements the computational
  // effort may be high, since the probability of getting a hole is
  // very low. 
  nvrtx = (nelemfat/2 > max_partgraph_vertices ?
	   max_partgraph_vertices : nelemfat);
  
  // Create adjacency table for partitioning with Metis. In the
  // adjacency graph the nodes are elements or group of elments of the
  // FEM mesh. Two nodes of the graph (elements of the mesh) have are
  // linked if they share a node.

  // adjncy:= xadj:= graph described in CSR format (as defined in
  // Metis documentation)
  int *xadj = new int[nvrtx+1],*adjncy;
  // vwgt:= weights for the vertices (elements) of the graph
  int *vwgt = new int[nvrtx];
  // el2vrtx:= maps elements to vertices when coalescing
  int *el2vrtx = new int[nelemfat];
  // vpart:= partitioning of vertices
  int *vpart = new int[nvrtx];
  vector<int> ngbrs_v;

  // Initialize 
  for (j=0; j<nelemfat; j++) el2vrtx[j] = -1;
  for (j=0; j<nvrtx; j++) vwgt[j] = 0;
  // Take `nvrtx' elements to coalesce the neighbor elements

  vrtx=0;
  while (!ngbrs.empty()) ngbrs.pop();

  if (nelemfat!=nvrtx) {
    // Take nvrtx `seeds' for growing the groups 
    while(1) {
      elem = irand(nelemfat);
      if (el2vrtx[elem]!=-1) continue;
      el2vrtx[elem]=vrtx++;
      vwgt[vrtx] += int(weight(elem)/weight_scale);
      ngbrs.push(elem);
      if (vrtx==nvrtx) break;
    }
    visited = nvrtx;
  } else {
    // Take a group for each element
    for (elem=0; elem<nelemfat; elem++) {
      el2vrtx[elem] = elem;
      find_elem(elem,nelemsetptr,nelemsets,mesh,elemset,locel);
      vwgt[elem] += int(weight(elem)/weight_scale);
    }      
    visited=nelemfat;
  }
  
  // Start coalescing neighbor elements until all elements are
  // assigned a vertex
  while (1) {
    if (ngbrs.empty()) {
      // either we have visited all the elements or there are
      // disconnected parts of the domain.

      // All elements have been visited
      if (visited==nelemfat) break;
      
      // Disconnected parts. Find next element not visited.
      for (k=0; k<nel; k++) {
	if (el2vrtx[k]<0) {
	  ngbrs.push(k);
	  // Put it in a random vertex
	  vrtx = irand(0,nvrtx)-1;
	  el2vrtx[k] = vrtx;
	  find_elem(k,nelemsetptr,nelemsets,mesh,elemset,locel);
	  vwgt[vrtx] += int(elemset->weight()/weight_scale);
	  break;
	}
      }
    }

    // Take the first element
    elem = ngbrs.front();
    ngbrs.pop();
    vrtx = el2vrtx[elem];
    assert(vrtx>=0);
    // Get its connectivities
    ngbrs_v.clear();
    set_ngbrs(elem,ngbrs_v);
    // Loop over all ngbrs of the fine vertex
    ne = ngbrs_v.end();
    for (n=ngbrs_v.begin(); n!=ne; n++) {
      int &elemk = *n;
      // Already visited
      if (el2vrtx[elemk]>=0) continue;
      // Assign to group
      el2vrtx[elemk]=vrtx;
      visited++;
      // Find elemset and add weight to vertex weight (for load
      // balancing) 
      vwgt[vrtx] += int(weight(elem)/weight_scale);
      // Put in queue
      ngbrs.push(elemk);
    }
  }

#if 1
  if (myrank==0) {
    for (elem=0; elem < nelemfat; elem++) {
      printf("el2vrtx[%d] = %d\n",elem,el2vrtx[elem]);
    }
  }
#endif

  // mark:= auxiliary vector that flags if an element has been marked
  // already as a linked node in the graph
#if 1
  int *mark = new int[nelemfat];
  for (int ielgj=0; ielgj<nelemfat; ielgj++) mark[ielgj]=-1;
#endif
        
  adjncy_v.resize(nvrtx);

  // Define graph descriptors 
  // Initializa graph ptrs
  for (elem=0; elem<nelemfat; elem++) {
    vrtxj = el2vrtx[elem];
    ngbrs_v.clear();
    set_ngbrs(elem,ngbrs_v);
    for (n=ngbrs_v.begin(); n!=ne; n++) {
      int &elemk = *n;
      vrtxjj = el2vrtx[ielgjj];
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

#if 1
  // print the graph
  for (int ielgj=0; ielgj<nvrtx; ielgj++) {
    printf("%d: ",ielgj);
    for (int jj=xadj[ielgj]; jj<xadj[ielgj+1]; jj++) 
      printf(" %d",adjncy[jj]);
    printf(", weight: %d\n",vwgt[ielgj]);
  }
#endif

  int options=0,edgecut,numflag=0,wgtflag=2;
  float *tpwgts_d = new float[npart];

#if 0
  for (j=0; j<nvsubdo; j++) 
    tpwgts_d[j] = tpwgts[subd2proc[j]/iisd_subpart]/float(iisd_subpart);
  if (size*iisd_subpart > 1) {
    if (myrank==0) {
      if (partflag==0) {
	if (myrank==0) printf("METIS partition - partflag = %d\n",partflag);
	METIS_WPartGraphKway(&nvrtx,xadj,adjncy,vwgt, 
			     NULL,&wgtflag,&numflag,&nvsubdo, 
			     tpwgts_d,&options,&edgecut,vpart);
      
      } else { // partflag=2
	assert(0);
      }
    }
  }  
#endif
}

